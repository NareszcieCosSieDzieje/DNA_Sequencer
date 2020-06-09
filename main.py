from sys import argv
import networkx as nx
import time
import matplotlib.pyplot as plt
import random


def check_overlap(seq: str, seq2: str):
    init_overlap = len(seq)
    for i in reversed(range(1, len(seq2))):
        if seq.endswith(seq2[:i]):
            return i
    return 0



def getSuccessor(graph, node):
    successors = graph.successors(node)
    sortedSuccessors = list(sorted(successors, key=lambda x: -graph[node][x]['weight']))
    mappedSuccessors = list(map(lambda x: (x, graph[node][x]['weight']), sortedSuccessors))
    # filteredSuccessors = list(filter(lambda x: x[1] == mappedSuccessors[0][1], mappedSuccessors))
    
    return mappedSuccessors


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

    
    start = time.time()
    # Create a graph, each node is a one k-gram
    content = []
    with open(argv[1]) as file:
        content = file.read().split('\n')[:-1]
    graph = nx.DiGraph()
    graph.add_nodes_from(content)
    
    edges = 0
    # Graph creation loop
    for seq in content:
            content2 = content[:]
            content2.remove(seq)
            for seq2 in content2:
                if ((overlap_count := check_overlap(seq, seq2)) > 0 ):
                    graph.add_edge(seq, seq2, weight = overlap_count)
                    edges += 1


    # Pruning the graph
    removedEdges = 0
    for n in graph.nodes():
        friends = graph.successors(n)
        maxWeight = 0
        for f in friends:
            currentWeight = graph[n][f]['weight']
            if currentWeight > maxWeight:
                maxWeight = currentWeight

        friends = graph.successors(n)
        edgesToRemove = []
        for f in friends:
            currentWeight = graph[n][f]['weight']
            if (currentWeight < maxWeight):
                removedEdges += 1
                edgesToRemove.append((n, f))
        graph.remove_edges_from(edgesToRemove)

    # print(edges, removedEdges, edges-removedEdges)

    eachNodeResults = []
    for n2 in graph.nodes():
        ourStack = []
        results = []

        path2 = [n2]
        iterWeight = 0
        iterEdges = 0
       
        #tuple = successor, weight
        futureSuccessors = getSuccessor(graph, n2)
        if len(futureSuccessors) < 1:
            break
        chosen = futureSuccessors[0]
        futureSuccessors.pop(0)
        for successor in futureSuccessors:
            ourStack.append({"node": successor, "weight": iterWeight, "edges": iterEdges, "path":path2})
        
        # pseudo Do While start
        while (len(list(graph.successors(chosen[0]))) > 0 and chosen[0] not in path2):
            # Add good successor to path
            path2.append(chosen[0])
            iterWeight += chosen[1]
            iterEdges += 1

            # Check successor, forks and add to stack
            futureSuccessors = getSuccessor(graph, chosen[0]) 
            if len(futureSuccessors) < 1:
                break
            chosen = futureSuccessors[0]
            futureSuccessors.pop(0)
            for successor in futureSuccessors:
                # print("Succ:", successor)
                ourStack.append({"node": successor, "weight": iterWeight, "edges": iterEdges, "path":path2})

        # Add to results array
        results.append([path2, iterWeight, iterEdges])

        while (ourStack):
            # Get state from stack and fill out things
            chosenState = ourStack[-1]
            ourStack.pop(-1)
            path2 = chosenState["path"]
            iterWeight = chosenState["weight"]
            iterEdges = chosenState["edges"]
            chosen = chosenState["node"]
            # print(ourStack)

            while (len(list(graph.successors(chosen[0]))) > 0 and chosen[0] not in path2):
                # Add good successor to path
                path2.append(chosen[0])
                iterWeight += chosen[1]
                iterEdges += 1

                # Check successor, forks and add to stack
                futureSuccessors = getSuccessor(graph, chosen[0]) 
                if len(futureSuccessors) < 1:
                    break
                chosen = futureSuccessors[0]
                futureSuccessors.pop(0)
                for successor in futureSuccessors:
                    ourStack.append({"node": successor, "weight": iterWeight, "edges": iterEdges, "path":path2})

            results.append([path2, iterWeight, iterEdges])
        # pseudo Do while end
    
        # Sum up results pick the best one go with the next iteration
        rpath = []
        rweight = 0
        redges = 0
        for r in results:
            if r[1] > rweight:
                rpath = r[0]
                rweight = r[1]
                redges = r[2]
        
        eachNodeResults.append([rpath, rweight, redges])

    # print(list(sorted(eachNodeResults, key=lambda x: x[1])))
    for j in reversed(list(sorted(eachNodeResults, key=lambda x: x[1]))):
        path = j[0]
        result = ""
        for i, node in enumerate(path):
            if i == 0:
                result = node
            else:
                overlapSlice = check_overlap(path[i-1], node)
                result += node[ overlapSlice :]
        print("Weight", "\t", "Edges", "\t", "Length", "\t", "Sequence")
        print(j[1], "\t", j[2], "\t", len(result), "\t", result)
        break

    stop = time.time()
    print("Execution time:", str((stop - start)*1000) + "ms\n")


    # Search for best path part
    # maxWeight = 0
    # maxResult = ""
    # maxEdgeCount = 0
    # maxResultLength = 0
    # path = []
    # # TODO: check what happends if you try dfs, maybe it will work???
    # for node in graph.nodes():
    #     path = [node]
    #     edgeCount = 0
    #     weight = 0
    #     # Get the successor with the best overlap
    #     succList = getSuccessor(graph, node)
    #     nextWeight = 0
    #     nextNode = succList[0][0]
    #     nextWeight = succList[0][1]

    #     # Starting from the successor add successors with the best overlap to the path
    #     # End when we have come through a cycle, each time increase weight and edgeCount
    #     lastBranch = None
    #     # a List of Lists containing []
    #     prevBranchStack = list()
    #     prevBranchWeight = 0
    #     # If there were more successors, for each successor append it to the list
    #     # together with weight value and edge count it had at that moment
    #     if len(succList) > 1:
    #         for i in succList[1:]:
    #             prevBranchStack.append({"node":i[0], "weight":i[1], "edgeCount":1})

    #     while (nextNode not in path and prevBranchStack == []):
    #         if nextNode in path:
    #             # to cofnij jestli jest prevBranch
    #             pass                
    #         path.append(nextNode)
    #         edgeCount += 1
    #         weight += nextWeight
    #         succList = getSuccessor(graph, nextNode)
    #         if len(succList) > 1:
    #             for i in succList[1:]:
    #                 prevBranchStack.append({"node":i[0], "weight":weight, "edgeCount":edgeCount})
    #         # if (succList[0][0] in path and len(succList) > 1):
    #         #     succList.pop(0)
           
    #         nextWeight = 0
    #         nextNode = succList[0][0]
    #         nextWeight = succList[0][1]
        

    #     print("nextNode Analysis: ", nextNode, list(graph.successors(path[-1])), [i in path for i in list(graph.successors(path[-1]))])
        
    #     # Walk through the path and concatenate k-mers to create a full sequence
    #     result = ""
    #     for i, node in enumerate(path):
    #         if i == 0:
    #             result = node
    #         else:
    #             overlapSlice = check_overlap(path[i-1], node)
    #             result += node[ overlapSlice :]
    #     print(weight, "\t", edgeCount, "\t", len(result), "\t", result)

    #     # The best solution will have most edges traversed and biggest weight
    #     if (edgeCount >= maxEdgeCount):
    #         if (weight > maxWeight):
    #             maxEdgeCount = edgeCount
    #             maxWeight = weight
    #             maxResult = result
    #             maxResultLength = len(maxResult)
    
    # # Print out the best result
    # print(maxResultLength, maxResult, maxEdgeCount)

   
    # red_edges = list(filter(lambda x: x[2] > 0, graph.edges(data='weight')))
    # #black_edges = list(filter(lambda x: x[2] <= 7, graph.edges(data='weight')))
    # path_edges = []
    # for i in range(len(path)-1):
    #     path_edges.append((path[i], path[i+1]))
    # blue_edges = path_edges
    # edge_labels = dict([((u, v,), d['weight']) for u, v, d in filter(lambda x: (x[0], x[1]) in path_edges, graph.edges(data=True))])

    # pos = nx.spring_layout(graph)
    # nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size = 500)
    # nx.draw_networkx_labels(graph, pos)
    # nx.draw_networkx_edges(graph, pos, edgelist=red_edges, edge_color='r', arrows=True)
    # nx.draw_networkx_edges(graph, pos, edgelist=blue_edges[0:4], edge_color='g', arrows=True)
    # nx.draw_networkx_edges(graph, pos, edgelist=blue_edges[4:], edge_color='b', arrows=True)
    # nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
    # plt.show()

    graph.clear()