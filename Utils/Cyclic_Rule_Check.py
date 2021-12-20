import sys

ChildrenVertices = dict()
root = None


def Get_Root():
    currs = []
    children_set = set()
    for curr, children in ChildrenVertices.items():
        currs.append(curr)
        for and_children in children:
            for child in and_children:
                children_set.add(child)
    roots = set()
    for curr in currs:
        if curr not in children_set:
            roots.add(curr)
    root_count = len(roots)
    if root_count != 1:
        raise Exception("multiple or no root(s), root count: {}".format(root_count))
    else:
        global root
        root = list(roots)[0]
        print('root: {}'.format(root))


def Cyclic_Rule_Check_Wrapper():
    result = Cyclic_Rule_Check(root, set([root]), [root])
    if result == []:
        print('no cycle')
    else:
        print('cycle found: {}'.format(result))


def Cyclic_Rule_Check(curr, nodes_set, nodes_list):
    if curr == -1:
        return []
    for and_children in ChildrenVertices[curr]:
        for child in and_children:
            # print(child)
            if child in nodes_set:
                print('[Caution]: cyclic rule found !!!')
                loop_start = nodes_list.index(child)
                return nodes_list[loop_start:]
            else:
                sub_nodes_set = nodes_set.copy()
                sub_nodes_set.add(child)
                sub_nodes_list = nodes_list.copy()
                sub_nodes_list.append(child)
                # print(sub_nodes_list)
                sub_return = Cyclic_Rule_Check(child, sub_nodes_set, sub_nodes_list)
                if sub_return:
                    return sub_return
    return []


def main():
    if len(sys.argv) != 2:
        print("usage: python3 Cyclic_Rule_Check.py [learned_tree file]")
        exit()

    global ChildrenVertices

    learned_tree = sys.argv[1]
    with open(learned_tree, "r") as tree_f:
        for line in tree_f:
            line = line.strip()
            tokens = line.split(',')
            children_count = int(tokens[0])
            curr_id = int(tokens[2])
            assert(children_count == (len(tokens) - 4) / 2)
            children_ids = []
            for j in range(children_count):
                index = 4 + 2 * j
                child_id = int(tokens[index])
                children_ids.append(child_id)
            if curr_id in ChildrenVertices:
                ChildrenVertices[curr_id].append(children_ids)
            else:
                ChildrenVertices[curr_id] = [children_ids]

    print(ChildrenVertices)

    Get_Root()

    Cyclic_Rule_Check_Wrapper()






if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Must use Python 3")
    main()
