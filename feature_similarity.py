from SPARQLWrapper import SPARQLWrapper, JSON
import os
import errno
import time
import numpy as np
from scipy.sparse import csr_matrix
import scipy
from scipy import spatial
from datetime import datetime
import pickle


class DomainMatrix:

    def __init__(self, endpoint, type_list, additional_filters, feature_strategy):

        print("{} Initialization of DomainMatrix...".format(datetime.now()))

        f_index_file = './tmp/f_index'
        mt_file = './tmp/sparse_matrix.npz'
        s_index_file = './tmp/s_index'

        if os.path.exists(mt_file):
            print("{} loading from files...".format(datetime.now()))
            self.mt = scipy.sparse.load_npz(mt_file)
            # load list of subject
            self.s_list = list()
            if os.path.exists(s_index_file):
                with open(s_index_file, 'rb') as f:
                    reader = f.readlines()
                    for l in reader:
                        self.s_list.append(l.strip())
            # based on test file, remove '
            self._rm_apostrophe()

            self.f_list = list()
            self.po = list()
            if os.path.exists(f_index_file):
                with open(f_index_file, 'rb') as f:
                    reader = f.readlines()
                    for i, l in enumerate(reader):
                        self.f_list.append(l.strip())
                        if '\t' in l or l.strip():
                            self.po.append(i)
            print("{} len of f_list: {}, len of po list: {}".format(datetime.now(), len(self.f_list), len(self.po)))
        else:
            print("creating new files...")
            sparql = SPARQLWrapper(endpoint)

            # prepare string of type list
            # for appending SPARQL
            if '<' not in type_list[0]:
                for i, t in enumerate(type_list):
                    t = '<' + t + '>'
                    type_list[i] = t
            str_type_list = ','.join(type_list)

            # prepare string of additional fitlers
            # for appending SPARQL
            str_filter_list = ' '.join(additional_filters)

            # ================== Subjects
            # ============================================

            # load list of subject
            self.s_list = list()
            if os.path.exists(s_index_file):
                with open(s_index_file, 'rb') as f:
                    reader = f.readlines()
                    for l in reader:
                        self.s_list.append(l.strip())
            else:   # create subject indices
                count_s_query = """ 
                                 PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                                 SELECT COUNT(DISTINCT ?s) as ?c
                                 WHERE { ?s rdf:type ?t .
                                 FILTER (?t IN (""" + str_type_list + """))}
                                 """

                sparql.setQuery(count_s_query)
                sparql.setReturnFormat(JSON)
                results = sparql.query().convert()
                for result in results["results"]["bindings"]:
                    count_s = int(result["c"]["value"])
                print("count of subjects: {}".format(count_s))

                # due to the limit of endpoint
                # we need several queries
                num_query = count_s // 10000 + 1
                for i in xrange(num_query):
                    q1 = """ 
                         PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                         SELECT distinct ?s
                         WHERE { ?s rdf:type ?t .
                         FILTER (?t IN (""" + str_type_list + """))}
                         LIMIT 10000 OFFSET """ + str(i * 10000)

                    sparql.setQuery(q1)
                    sparql.setReturnFormat(JSON)
                    results = sparql.query().convert()
                    try:
                        os.makedirs(os.path.dirname(s_index_file))
                    except OSError as exc:  # Guard against race condition
                        if exc.errno != errno.EEXIST:
                            raise
                    for result in results["results"]["bindings"]:
                        self.s_list.append(result["s"]["value"].encode("UTF-8"))
                    print("{} finished {}%\r".format(datetime.now(), 100.0 * i / num_query))

                with open(s_index_file, 'a') as f:
                    for sub in self.s_list:
                        f.write(sub + '\n')
            # self.s_list = list(set(self.s_list))
            print('total subjects: {}'.format(len(self.s_list)))

            print("")

            # ================== Build Matrix & Features
            # ============================================
            if os.path.exists('./tmp/col.pkl'):
                col = pickle.load('./tmp/col.pkl')
                row = pickle.load('./tmp/row.pkl')
            else:
                row, col = list(), list()

            # load list of features
            # self.f_list = list()
            # in case a public endpoint not responding
            # and restart to continue previous retrieval

            if os.path.exists(f_index_file):
                with open(f_index_file, 'rb') as f:
                    reader = f.readlines()
                    for l in reader:
                        self.f_list.append(l.strip())
            else:
                self.f_list = list()

            continue_offset = 0

            # ======= incoming node-predicate as features
            # ===================================================
            if 'subject-predicate' in feature_strategy:

                count_sp_query = """ 
                                 PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                                 SELECT (COUNT(*) AS ?c) {
                                    SELECT ?s ?i ?j 
                                    WHERE { ?s rdf:type ?t . ?i ?j ?s . 
                                    FILTER (?t IN (""" + str_type_list + """))
                                    """ + str_filter_list + """
                                }}
                                """

                sparql.setQuery(count_sp_query)
                sparql.setReturnFormat(JSON)
                results = sparql.query().convert()
                for result in results["results"]["bindings"]:
                    count_sp = int(result["c"]["value"])
                print("count of incoming triples: {}".format(count_sp))

                # due to the limit of endpoint
                # we need several queries
                num_query = count_sp // 10000 + 1
                # num_query = 1
                for nq in xrange(num_query):
                    if nq >= continue_offset:
                        q1 = """ 
                             PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                             SELECT ?s ?i ?j 
                             WHERE { ?s rdf:type ?t . ?i ?j ?s .
                             FILTER (?t IN (""" + str_type_list + """)) 
                             FILTER (!isLiteral(?z)) """ + str_filter_list + """
                             }
                             LIMIT 10000 OFFSET """ + str(i * 10000)

                        sparql.setQuery(q1)
                        sparql.setReturnFormat(JSON)
                        results = sparql.query().convert()
                        try:
                            os.makedirs(os.path.dirname(f_index_file))
                        except OSError as exc:  # Guard against race condition
                            if exc.errno != errno.EEXIST:
                                raise
                        for result in results["results"]["bindings"]:
                            s = result["s"]["value"].encode("UTF-8")
                            i = result["i"]["value"].encode("UTF-8")
                            j = result["j"]["value"].encode("UTF-8")

                            # predicate+object as a feature
                            f = i + '\t' + j
                            row.append(self.s_list.index(s))
                            if f in self.f_list:
                                col.append(self.f_list.index(f))
                            else:
                                self.f_list.append(f)
                                col.append(len(self.f_list) - 1)

                        print("{} finished {}%, OFFSET={}\r".format(datetime.now(), 100.0 * nq / num_query, i * 10000))
                        time.sleep(30)

                        # write row/col to pickle
                        with open('./tmp/row.pkl', 'wb') as pf:
                            pickle.dump(row, pf)
                        with open('./tmp/col.pkl', 'wb') as pf:
                            pickle.dump(col, pf)

                        # write features to file
                        with open(f_index_file, 'wb') as f:
                            for ft in self.f_list:
                                f.write(ft + '\n')

            # ============= outgoing objects
            # ===============================================
            count_f_query = """ 
                             PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                             SELECT (COUNT(*) AS ?c) {
                                SELECT ?s ?y ?z 
                                WHERE { ?s rdf:type ?t . ?s ?y ?z . 
                                FILTER (?t IN (""" + str_type_list + """))
                                FILTER (!isLiteral(?z)) """ + str_filter_list + """
                            }}
                            """

            sparql.setQuery(count_f_query)
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()
            for result in results["results"]["bindings"]:
                count_f = int(result["c"]["value"])
            print("count of outgoing triples: {}".format(count_f))

            # due to the limit of endpoint
            # we need several queries
            num_query = count_f // 10000 + 1
            # num_query = 1
            for i in xrange(num_query):
                if i >= continue_offset:
                    q1 = """ 
                         PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                         SELECT ?s ?y ?z 
                         WHERE { ?s rdf:type ?t . ?s ?y ?z .
                         FILTER (?t IN (""" + str_type_list + """)) 
                         FILTER (!isLiteral(?z)) """ + str_filter_list + """
                         }
                         LIMIT 10000 OFFSET """ + str(i * 10000)

                    sparql.setQuery(q1)
                    sparql.setReturnFormat(JSON)
                    results = sparql.query().convert()
                    try:
                        os.makedirs(os.path.dirname(f_index_file))
                    except OSError as exc:  # Guard against race condition
                        if exc.errno != errno.EEXIST:
                            raise
                    for result in results["results"]["bindings"]:
                        s = result["s"]["value"].encode("UTF-8")
                        y = result["y"]["value"].encode("UTF-8")
                        z = result["z"]["value"].encode("UTF-8")

                        # predicate+object as a feature
                        f = y + '\t' + z
                        row.append(self.s_list.index(s))
                        if f in self.f_list:
                            col.append(self.f_list.index(f))
                        else:
                            self.f_list.append(f)
                            col.append(len(self.f_list)-1)

                        if 'object' in feature_strategy:
                            # object as a feature
                            row.append(self.s_list.index(s))
                            if z in self.f_list:
                                col.append(self.f_list.index(z))
                            else:
                                self.f_list.append(z)
                                col.append(len(self.f_list)-1)

                            # subject as a feature
                            row.append(self.s_list.index(s))
                            if s in self.f_list:
                                col.append(self.f_list.index(s))
                            else:
                                self.f_list.append(s)
                                col.append(len(self.f_list)-1)

                    print("{} finished {}%, OFFSET={}\r".format(datetime.now(), 100.0*i/num_query, i*10000))
                    time.sleep(30)

                    # write row/col to pickle
                    with open('./tmp/row.pkl', 'wb') as pf:
                        pickle.dump(row, pf)
                    with open('./tmp/col.pkl', 'wb') as pf:
                        pickle.dump(col, pf)

                    # write features to file
                    with open(f_index_file, 'wb') as f:
                        for ft in self.f_list:
                            f.write(ft + '\n')

            print('current total features: {}'.format(len(self.f_list)))

            # construct matrix
            col = np.asarray(col)
            row = np.asarray(row)
            data = np.ones(shape=np.shape(col), dtype=np.float)
            self.mt = csr_matrix((data, (row, col)), shape=(len(self.s_list), len(self.f_list)))
            scipy.sparse.save_npz(mt_file, self.mt)

            # based on test file, remove ' in subjects
            self._rm_apostrophe()

    def _rm_apostrophe(self):
        for idx, s in enumerate(self.s_list):
            self.s_list[idx] = s.replace('\'', '')

    def rec(self, test_file, filter_p_list=None):
        """Provide recommendation

        Args:
            test_file: tes file
            filter_p_list: properties to be filtered
            in the f_index while providing recommendations
        """
        # get direct relationship features
        # get object of filtering properties
        dir_f_ind, filter_p_list_o, filter_p_list_ind = \
            list(), list(), list()
        with open('./tmp/f_index', 'r') as f:
            reader = f.readlines()
            for ind, val in enumerate(reader):
                if '\t' not in val:
                    dir_f_ind.append(ind)
                else:
                    if filter_p_list is not None:
                        splitted = val.strip().split('\t')
                        if splitted[0] in filter_p_list:
                            filter_p_list_o.append(splitted[1])
                            filter_p_list_ind.append(ind)

        print("{} direct node features: {}"
              .format(datetime.now(), len(dir_f_ind)))
        
        if filter_p_list is not None:
            filter_p_list_o = set(filter_p_list_o)
            print("{} filtering property node features: {}"
              .format(datetime.now(), len(filter_p_list_o)))

            with open('./tmp/f_index', 'r') as f:
                reader = f.readlines()
                for ind, val in enumerate(reader):
                    if '\t' not in val and val.strip() in filter_p_list_o:
                        filter_p_list_ind.append(ind)
            print("{} filtering features in total {}"
                  .format(datetime.now(), len(filter_p_list_ind)))
                    
        # np_mt = self.mt.toarray()   # Memory error
        # print("{} converted sparse to numpy array".format(datetime.now()))
        # col_sum = np.sum(np_mt, axis=0)
        # log_idf_np_mt = np.log(col_sum/(np_mt+1.0))
        # binary format of mt 0, 1
        self.mt = self.mt.sign()
        col_sum = self.mt.sum(axis=0)
        print("{} type(col_sum): {}, shape: {}"
              .format(datetime.now(), type(col_sum), np.shape(col_sum)))
        col_count = self.mt.getnnz(axis=1)
        print("type(col_count): {}, shape: {}"
              .format(type(col_count), np.shape(col_count)))
        print("{} triple has features in std: {}, mean:{}"
              .format(datetime.now(), np.std(col_count), np.mean(col_count)))
        mean_col_count = np.mean(col_count)
        idf = np.log(len(self.s_list)/(col_sum+1.0))
        idf = np.squeeze(np.asarray(idf))
        # discount direct node features,
        # 0.0 means no direct nodes
        discount = 1.0
        idf[dir_f_ind] = idf[dir_f_ind] * discount
        print("{} idf shape: {}".format(datetime.now(), np.shape(idf))) #3867139
        # log_mt = self.mt.log1p()

        # use mask to select
        # property related ones
        if filter_p_list is not None:
            mask = np.zeros(len(idf))
            mask[filter_p_list_ind] = 1.0
            idf = idf * mask

        # generate candidates, which are all
        # items in the test file
        candidates = list()
        with open(test_file, 'rb') as f:
            reader = f.readlines()
            for l in reader:
                splits = l.strip().split('\t')
                item = splits[0]
                items = splits[1].strip().split(' ')
                if item not in self.s_list:
                    print("target item {} not in the subject list".format(item))
                candidates.append(item)
                for i in items:
                    if i not in self.s_list:
                        print("items {} not in the subject list".format(i))
                    else:
                        candidates.append(i)
        candidates = list(set(candidates))
        print("{} candidate items: {}".format(datetime.now(), len(candidates)))

        np_c = list()
        c_indices = list()
        for c in candidates:
            c_ind = self.s_list.index(c)
            c_indices.append(c_ind)
            b = self.mt.getrow(c_ind).nonzero()[1]
            # print b
            np_c.append(b)
        print("{} candidate np_c len: {}".format(datetime.now(), len(np_c)))

        with open(test_file, 'rb') as f:
            reader = f.readlines()
            for idx, l in enumerate(reader):
                candidates_scores = dict()
                if idx % 50 == 0:
                    print("{} {}-th line".format(datetime.now(), idx))
                splits = l.strip().split('\t')
                item = splits[0]
                rec_str = item + '\t'
                item_ind = self.s_list.index(item)
                try:
                    # optimized version
                    a = self.mt.getrow(item_ind).nonzero()[1]
                    # beside IDF part
                    t1 = (1.2 + 1.0) / (1.0 + 1.2 * (1.0 - 0.75 + 0.75 * col_count[c_indices] / mean_col_count)) + 1.0
                    # print("{} t1 shape: {}".format(datetime.now(), np.shape(t1)))
                    # idf part
                    for c_ind, c in enumerate(candidates):
                        overlap = [x for x in a if x in np_c[c_ind]]
                        t2 = idf[overlap]
                        # print("{} t2 shape: {}".format(datetime.now(), np.shape(t2)))
                        s = np.sum(t1[c_ind] * t2)
                        candidates_scores[c] = s

                    # initial version
                    # for c in candidates:
                    #     c_ind = self.s_list.index(c)
                    #     a = self.mt.getrow(item_ind).toarray().reshape((-1,))
                    #     b = self.mt.getrow(c_ind).toarray().reshape((-1,))
                    #     # only po
                    #     # a = self.mt.getrow(item_ind).toarray().reshape((-1,))[self.po]
                    #     # b = self.mt.getrow(c_ind).toarray().reshape((-1,))[self.po]
                    #
                    #     # BM25+, po-list
                    #     # BM25+ parameters
                    #     # k1 = 1.2
                    #     # b = 0.75
                    #     # delta = 1.0
                    #     # overlap = np.nonzero(a * b)
                    #     # s = np.sum(idf[overlap] * (
                    #     #         (1.2 + 1.0) / (1.0 + 1.2 * (1.0 - 0.75 + 0.75 * col_count[c_ind] / mean_col_count))
                    #     #         + 1.0
                    #     #     ))
                    #
                    #     # TFIDF, po-list, cosine
                    #     # s = 1.0 - spatial.distance.cosine((a * idf)[self.po], (b * idf)[self.po])
                    #
                    #     # TFIDF, cosine
                    #     # s = 1.0 - spatial.distance.cosine((a * idf), (b * idf))
                    #
                    #     # TFIDF, smoothed cosine
                    #     # overlap = len(np.nonzero(a * b))
                    #     # smoothing = (overlap / (overlap + 10.0))
                    #     # s = smoothing * (1.0 - spatial.distance.cosine((a * idf), (b * idf)))
                    #
                    #     # PICSS
                    #     s = 1.0 - spatial.distance.jaccard(a, b, idf)

                        # non_zero_a = np.nonzero(a)
                        # non_zero_b = np.nonzero(b)
                        # intersect = np.intersect1d(non_zero_a, non_zero_b)
                        # union = np.union1d(non_zero_a, non_zero_b)
                        # s = np.sum(idf[intersect]) / np.sum(idf[union])
                        candidates_scores[c] = s
                except KeyError:
                    print("item {} keyerror".format(item))

                # sort and write
                for key, value in reversed(sorted(candidates_scores.iteritems(), key=lambda (k, v): (v, k))):
                    rec_str += (key + ':' + str(value) + ' ')
                rec_str += '\n'
                with open('rec_FSSM_discount1.0_top15', 'a') as f:
                    f.write(rec_str)
                # print("{} sort scores and write to the file...".format(datetime.now()))


if __name__ == '__main__':

    endpoint = 'http://dbpedia.org/sparql'
    # endpoint = 'http://140.203.154.138:8018/sparql'
    type_list = ['http://dbpedia.org/ontology/Book']
    feature_strategy = ['predicate-object', 'subject-predicate']

    # filters
    # filter_z = "FILTER STRSTARTS(STR(?z), 'http://dbpedia.org/resource')"
    # filter_y = "FILTER (?p IN (<http://purl.org/dc/terms/subject>, " \
    #            "<http://dbpedia.org/ontology/author>, " \
    #            "<http://dbpedia.org/ontology/publisher>, " \
    #            "<http://dbpedia.org/ontology/literaryGenre>, " \
    #            "<http://dbpedia.org/ontology/mediaType>, " \
    #            "<http://dbpedia.org/ontology/subsequentWork>, " \
    #            "<http://dbpedia.org/ontology/previousWork>, " \
    #            "<http://dbpedia.org/ontology/country>, " \
    #            "<http://dbpedia.org/ontology/series>, " \
    #            "<http://dbpedia.org/ontology/nonFictionSubject>, " \
    #            "<http://dbpedia.org/ontology/coverArtist>, " \
    #            "<http://dbpedia.org/ontology/illustrator>, " \
    #            "<http://dbpedia.org/ontology/genre>, " \
    #            "<http://dbpedia.org/ontology/translator>, " \
    #            "<http://dbpedia.org/ontology/recordLabel>))"
    filter_y = ["http://purl.org/dc/terms/subject",
               "http://dbpedia.org/ontology/genre",
               "http://dbpedia.org/ontology/associatedBand",
               "http://dbpedia.org/ontology/associatedMusicalArtist",
               "http://dbpedia.org/ontology/instrument",
               "http://dbpedia.org/ontology/recordLabel",
               "http://dbpedia.org/ontology/occupation",
               "http://dbpedia.org/ontology/hometown",
               "http://dbpedia.org/ontology/bandMember",
               "http://dbpedia.org/ontology/formerBandMember",
               "http://dbpedia.org/ontology/currentMember",
               "http://dbpedia.org/ontology/influencedBy",
               "http://dbpedia.org/ontology/pastMember",
               "http://dbpedia.org/ontology/associatedAct",
               "http://dbpedia.org/ontology/influenced"]
    additional_filters = list()
    # additional_filters.append(filter_z)
    # additional_filters.append(filter_y)

    DM = DomainMatrix(endpoint, type_list, additional_filters, feature_strategy)

    domain = 'music'
    test_file = '/Users/parklize/Documents/Workspaces/Java/EntityEmb/data/' + domain + '/test'
    # test_file = '/home/guapia/Program/EntityEmb/bin/data/' + domain + '/test_itemCountBT10_filtered'
    DM.rec(test_file)