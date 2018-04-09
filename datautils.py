

def filtered_test(test_file, s_index):

    s_list = list()
    with open(s_index, 'rb') as f:
        reader = f.readlines()
        for l in reader:
            s_list.append(l.strip())
    for idx, s in enumerate(s_list):
        s_list[idx] = s.replace('\'', '')

    with open(test_file, 'rb') as f:
        reader = f.readlines()
        for l in reader:
            splitted = l.strip().split('\t')
            target_item = splitted[0]
            items = splitted[1].split(' ')
            if target_item in s_list:
                wl = target_item + '\t'
                for it in items:
                    if it in s_list:
                        wl += (it + ' ')
                with open(test_file+'_filtered', 'a') as fw:
                    fw.write(wl.strip()+'\n')
            else:
                print("{} not in s_index".format(target_item))


def stats(test_file):
    rating_list = list()
    with open(test_file) as f:
        reader = f.readlines()
        for l in reader:
            splitted = l.strip().split('\t')
            target_item = splitted[0]
            items = splitted[1].split(' ')
            rating_list.append(len(items))
    print(float(sum(rating_list))/len(rating_list))
    print(sum(rating_list))


def modify_to_single_test(train, test):

    train_dict = dict()
    with open(train, 'rb') as f:
        reader = f.readlines()
        for line in reader:
            split = line.split('\t')
            train_dict[split[0].encode('UTF-8')] = split[1].strip().encode('UTF-8')

    with open(test, 'rb') as f:
        reader = f.readlines()
        for line in reader:
            split = line.split('\t')
            items = split[1].strip().split(' ')
            given_item = train_dict[split[0]]
            rec = given_item + '\t'
            for item in items:
                if item != given_item:
                    rec += (item.encode('UTF-8') + ' ')
            with open('facebook_500_test', 'a') as f:
                f.write(rec+'\n')


if __name__ == '__main__':
    # domain = 'book'
    # test_file = '/Users/parklize/Documents/Workspaces/Java/EntityEmb/data/' + domain + '/test_itemCountBT10'
    # s_index = './tmp/s_index'
    # filtered_test(test_file, s_index)
    # stats(test_file+'_filtered')

    # modify facebook train/test to a single test file
    # train = '/Users/parklize/Documents/Workspaces/Python2/Thesis/dataset/C5_FacebookLikesMusicDomain_500_train.tsv'
    # test = '/Users/parklize/Documents/Workspaces/Python2/Thesis/dataset/C5_FacebookLikesMusicDomain_500_test.tsv'
    # modify_to_single_test(train, test)
    test_file = 'dataset/test'
    s_index = './tmp/s_index'
    filtered_test(test_file, s_index)