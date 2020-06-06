from sys import argv
import networkx as nx
import matplotlib.pyplot as plt
import random


def check_overlap(seq: str, seq2: str):
    init_overlap = len(seq)
    for i in range(len(seq2)):
        if seq.endswith(seq2[:init_overlap]):
            return init_overlap
        init_overlap -= 1
    
    return init_overlap


def getSuccessor(graph, node):
    successors = graph.successors(node)
    sortedSuccessors = list(sorted(successors, key=lambda x: -graph[node][x]['weight']))
    mappedSuccessors = list(map(lambda x: (x, graph[node][x]['weight']), sortedSuccessors))
    filteredSuccessors = list(filter(lambda x: x[1] == mappedSuccessors[0][1], mappedSuccessors))
    
    return filteredSuccessors

def getSuccessorV2(graph, node):
    successors = graph.successors(node)
    bestWeight = 0
    bestSuccesor = ""
    bestSuccesor = random.choice(list(successors))
    bestWeight = graph[node][bestSuccesor]['weight']
    return bestSuccesor, bestWeight


def getSuccessorV3(graph, node):
    successors = graph.successors(node)
    bestSuccesor = ""
    bestSuccessors = dict()

    for successor in successors:
        currentWeight = graph[node][successor]['weight']

        if (len(bestSuccessors.keys()) < 2):
            bestSuccessors[successor] = currentWeight

    bestSuccesor = random.choice(list(bestSuccessors.keys()))
    bestWeight = bestSuccessors[bestSuccesor]

    return bestSuccesor, bestWeight


def getSuccessorV4(graph, node, butno = []):
    successors = graph.successors(node)
    sortedSuccessors = sorted(successors, key=lambda x: -graph[node][x]['weight'])
    sortedSuccessors = list(filter(lambda x: x not in butno, sortedSuccessors))
    if (len(sortedSuccessors) == 0):
        return "ERROR", 0

    return sortedSuccessors[0], graph[node][sortedSuccessors[0]]['weight']


if __name__ == "__main__":
    if (len(argv) < 2):
        print("If you want to recreate a sequence based on a spectrum of nucleotides:")
        print("Usage: " + argv[0] + " [filename]")
        print("If you want to create a spectrum from sequence and randomly delete some of the reads:")
        print("Usage: " + argv[0] + " [filename] [len_of_read] [delete_count]")
        exit(1)

    # Create a sequence file for testing
    if (len(argv) == 4):
        spectrum = []
        sequence = ""
        with open(argv[1]) as file:
            sequence = file.read()
        arglen = int(argv[2])
        deletecount = int(argv[3])
        for i in range(len(sequence) - arglen + 1):
            spectrum.append(sequence[i:i+arglen])

        for i in range(deletecount):
            spectrum.remove(random.choice(spectrum))
        random.shuffle(spectrum)
        for i in spectrum:
            print(i.upper())
        exit(0)
    
    # Create a graph, each node is a one k-gram
    content = []
    with open(argv[1]) as file:
        content = file.read().split('\n')[:-1]
    graph = nx.DiGraph()
    graph.add_nodes_from(content)
    
    
    # Graph creation loop
    gCopy = graph.copy()
    limit = 1
    flag = False
    while (True):
        # Create edges if overlap is bigger than limit
        edges = []
        for seq in content:
            content2 = content[:]
            content2.remove(seq)
            for seq2 in content2:
                if ((overlap_count := check_overlap(seq, seq2)) > limit):
                    edges.append([seq, seq2])
                    graph.add_edge(seq, seq2, weight = overlap_count)
        
        # Check if there is a node without a edge
        for n in graph.nodes():
            # If yes restore graph from copy and break the loop
            if len(list(graph.neighbors(n))) < 1:
                graph = gCopy
                flag = True
                break
        
        # Break the loop
        if (flag):
            break

        # If there isn't such a node, increase limit and save current graph then clear edges from it
        limit += 1
        gCopy = graph.copy()
        graph.remove_edges_from(edges)

    
    # Search for best path part
    maxWeight = 0
    maxResult = ""
    maxEdgeCount = 0
    maxResultLength = 0
    path = []
    # TODO: check what happends if you try dfs, maybe it will work???
    for node in graph.nodes():
        print(nx.dfs_successors(graph, node))
        path = [node]
        edgeCount = 0
        weight = 0
        # Get the successor with the best overlap
        succList = getSuccessor(graph, node)
        nextWeight = 0
        # if len(succList ) > 1:
        #     for succ in succList:
        #         succc = getSuccessor(graph, succ[0])
        #         if succc[0][1] > nextWeight:
        #             nextWeight = succ[1]
        #             nextNode = succ[0]
        # else:
        nextNode = succList[0][0]
        nextWeight = succList[0][1]
        # print(succList)
        # nextNode = succList[0][0]
        # nextWeight = succList[0][1]
        # print(path, nextNode)

        # Starting from the successor add successors with the best overlap to the path
        # End when we have come through a cycle, each time increase weight and edgeCount
        while (nextNode not in path):
            # print("in while")
            path.append(nextNode)
            # print(nextNode)
            edgeCount += 1
            weight += nextWeight
            # succList = getSuccessor(graph, nextNode)
            succList = getSuccessor(graph, nextNode)
            nextWeight = 0
            # if len(succList ) > 1:
            #     for succ in succList:
            #         succc = getSuccessor(graph, succ[0])
            #         if succc[0][1] > nextWeight:
            #             nextWeight = succ[1]
            #             nextNode = succ[0]
            # else:
            nextNode = succList[0][0]
            nextWeight = succList[0][1]
            # bestNextNode = succList[0][0]
            # nextWeight = succList[0][1]
            # nextNode = bestNextNode
        
        # Walk through the path and concatenate k-mers to create a full sequence
        result = ""
        for i, node in enumerate(path):
            if i == 0:
                result = node
                #print("Begin =", node)
            else:
                overlapSlice = check_overlap(path[i-1], node)
                result += node[ overlapSlice :]
                #print("Do dodania = {z} ,do uciecia = {x}, po ucieciu = {y}, wynik = {q}".format(x=node[0:overlapSlice], y=node[overlapSlice:], z=path[i-1], q=result))
        print(weight, "\t", edgeCount, "\t", len(result), "\t", result)
        #print(len(result))

        # The best solution will have most edges traversed and biggest weight
        if (edgeCount >= maxEdgeCount):
            if (weight > maxWeight):
                maxEdgeCount = edgeCount
                maxWeight = weight
                maxResult = result
                maxResultLength = len(maxResult)
    
    # Print out the best result
    print(maxResultLength, maxResult, maxEdgeCount)

   
    red_edges = list(filter(lambda x: x[2] > 8, graph.edges(data='weight')))
    black_edges = list(filter(lambda x: x[2] <= 7, graph.edges(data='weight')))
    path_edges = []
    for i in range(len(path)-1):
        path_edges.append((path[i], path[i+1]))
    blue_edges = path_edges
    edge_labels = dict([((u, v,), d['weight']) for u, v, d in filter(lambda x: (x[0], x[1]) in path_edges, graph.edges(data=True))])

    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size = 500)
    nx.draw_networkx_labels(graph, pos)
    nx.draw_networkx_edges(graph, pos, edgelist=red_edges, edge_color='r', arrows=True)
    nx.draw_networkx_edges(graph, pos, edgelist=blue_edges[0:1], edge_color='g', arrows=True)
    nx.draw_networkx_edges(graph, pos, edgelist=blue_edges[1:], edge_color='b', arrows=True)
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
    plt.show()

    graph.clear()