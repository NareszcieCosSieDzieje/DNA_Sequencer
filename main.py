from sys import argv
from sys import stderr
import networkx as nx
import time
import matplotlib.pyplot as plt
import random
from statistics import median
import numpy as np
import re


def check_overlap(seq: str, seq2: str):
    init_overlap = len(seq)
    for i in reversed(range(1, len(seq2))):
        if seq.endswith(seq2[:i]):
            return i
    return 0


def concatenatePath(path: list):
    result = ""
    for i, el in enumerate(path):
        if i == 0:
            result = el
        else:
            overlap = check_overlap(path[i-1], el)
            result += el[overlap:]
    return result


def trimGraph(graph: nx.DiGraph, aboveTreshPercentage: float=0.8, percentileLower: int=50, percentileUpper: int=90):
    """Prune the graph to reduce number of edges
    
    Keyword arguments:
        aboveTreshPercentage -- Percent of average weights that won't be removed (default 0.8)
        percentileLower -- Above this value edge may remain with aboveThreshPercentage chance (default 50)
        percentileUpper -- Above this value edge is sure to remain (default 90)
    """
    removedEdges = 0
    # This loop has O(n*(n+n)) = O(n^2)
    for n in graph.nodes():
        friends = graph.successors(n)

        weights = []
        for f in friends:
            currentWeight = graph[n][f]['weight']
            weights.append(currentWeight)
        
        limitWeightUpper = np.percentile(weights, percentileUpper)
        limitWeightLower = np.percentile(weights, percentileLower)

        # Create a list of edges to remove
        friends = graph.successors(n)
        edgesToRemove = []
        initialFriends = len(list(friends))
        friends = graph.successors(n)
        for f in friends:
            currentWeight = graph[n][f]['weight']
            # If weight is below median remove it
            if (currentWeight <= limitWeightLower and initialFriends > 1):
                initialFriends -= 1
                removedEdges += 1
                edgesToRemove.append((n, f))
            # Else if edge is above median but below sure threshold
            # it has aboveThreshPercentage to not be removed
            elif (currentWeight <= limitWeightUpper and initialFriends > 1):
                chance = random.uniform(0, 1)
                if (chance < aboveTreshPercentage):
                    initialFriends -= 1
                    removedEdges += 1
                    edgesToRemove.append((n, f))
            
        # Execute order 66
        graph.remove_edges_from(edgesToRemove)
    
    return removedEdges


def heuristicSucc(graph: nx.DiGraph, optimalValue: int, maxLength: int, looseness: int = 0):

    bestPath = []
    bestLen = 0
    for n in graph.nodes():
        path = [n]
        bestSucc = n
        pathLen = len(path[0])
        while (pathLen < maxLength and len(path) != optimalValue):
            successors = getSuccessor(graph, bestSucc)
            if len(successors) == 0:
                break

            # Find max weight
            maxWeight = 0
            for succW in successors:
                currWeight = graph[bestSucc][succW[0]]['weight']
                if currWeight > maxWeight:
                    maxWeight = currWeight
            
            # Get successors close to max weight
            succs = []
            for succ in successors:
                currWeight = graph[bestSucc][succ[0]]['weight']
                if (currWeight >= maxWeight - looseness):
                    succs.append(succ)
            if len(succs) == 0:
                break

            bestSucc = ""
            maxUniqueSuccs = 0
            for succ in succs:
                uniqueSuccs = 0
                successorSuccessors = getSuccessor(graph, succ[0])
                # print(successorSuccessors)
                if len(successorSuccessors) == 0:
                    continue

                for succSuccs in successorSuccessors:
                    if succSuccs[0] not in path:
                        uniqueSuccs += 1

                if uniqueSuccs >= maxUniqueSuccs:
                    bestSucc = succ[0]
                    maxUniqueSuccs = uniqueSuccs

            pathLen += len(path[0]) - check_overlap(path[-1], bestSucc)
            if pathLen > maxLength:
                pathLen -= len(path[0]) - check_overlap(path[-1], bestSucc)
                break
            path.append(bestSucc)
            # print("Unique succs:{0}".format(maxUniqueSuccs))
        # print("\nNode:{0}\nPath:{1}\npathLen:{2}\nlen(path):{3}".format(n, path, pathLen, len(path)))
        if len(set(path)) > bestLen:
            bestPath = [*path]
            bestLen = len(set(path))

    return bestPath, bestLen, optimalValue - bestLen



def heuristicSucc2(graph: nx.DiGraph, optimalValue: int, maxLength: int):

    bestPath = []
    bestLen = 0
    for n in graph.nodes():
        path = [n]
        bestSucc = n
        pathLen = len(path[0])
        while (pathLen < maxLength and len(path) != optimalValue):
            successors = getSuccessor(graph, bestSucc)
            if len(successors) == 0:
                break

           
            # bestSucc = ""
            
            bestUniqSuccs = dict()
            for succ in successors:
                uniqueSuccs = 0
                successorSuccessors = getSuccessor(graph, succ[0])
                # print(successorSuccessors)
                if len(successorSuccessors) == 0:
                    continue
              
                for succSuccs in successorSuccessors:
                    if succSuccs[0] not in path:
                        uniqueSuccs += 1

                bestUniqSuccs[succ[0]] = uniqueSuccs
            bestUniqSuccsSorted = [elem[0] for elem in sorted(bestUniqSuccs.items(), reverse=True, key=lambda x: x[1])][:3]

            if len(bestUniqSuccsSorted) == 0:
                break
           
            maxWeight = 0
            nextBestSucc = ""
            for succW in bestUniqSuccsSorted:
                currWeight = graph[bestSucc][succW]['weight']
                if currWeight > maxWeight:
                    maxWeight = currWeight
                    nextBestSucc = succW
                    
            bestSucc = nextBestSucc
            

            pathLen += len(path[0]) - check_overlap(path[-1], bestSucc)
            if pathLen > maxLength:
                pathLen -= len(path[0]) - check_overlap(path[-1], bestSucc)
                break
            path.append(bestSucc)
            
            # print("Unique succs:{0}".format(maxUniqueSuccs))
        # print("\nNode:{0}\nPath:{1}\npathLen:{2}\nlen(path):{3}".format(n, path, pathLen, len(path)))
        if len(set(path)) > bestLen:
            bestPath = [*path]
            bestLen = len(set(path))

    return bestPath, bestLen, optimalValue - bestLen


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


def findBestSequence():
    pass


if __name__ == "__main__":

    random.seed(time.time())
    #random.seed(5)
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
        for i in range(len(sequence) - arglen + 1):
            spectrum.append(sequence[i:i+arglen])

        if "+" in argv[3]: # Add random
            '''
            modifycount = int(argv[3][1:])
            added = 0
            while (added != modifycount):
                random_nuc = ""
                for i in range(arglen):
                    random_nuc += random.choice(["A","C","T","G"])
                if random_nuc not in spectrum:
                    spectrum.append(random_nuc)
                    added += 1
            '''
            modifycount = int(argv[3][1:])
            added = 0
            while (added != modifycount):
                random_nuc = ""
                random_nuc = list(random.choice(spectrum))
                randomElements = []
                for index, x in enumerate(random_nuc):
                    if random.randint(0, 1):
                        y = x
                        while(y == x):
                            y = random.choice(['A', 'C', 'T', 'G'])
                        random_nuc[index] = y
                random_nuc = "".join(random_nuc)
                if random_nuc not in spectrum:
                    spectrum.append(random_nuc)
                    added += 1
        else: # Delete random
            modifycount = int(argv[3])

            for i in range(modifycount):
                spectrum.remove(random.choice(spectrum))

        random.shuffle(spectrum)
        print(len(set(spectrum)), "==", len(spectrum), file=stderr)
        for i in spectrum:
            print(i.upper())
        exit(0)


    

    txt = argv[1]
    pattern = re.compile(".*\.(?P<numberix>[0-9]+)[+-](?P<extraSize>[0-9]+).*") 
    match = re.findall(pattern, txt)
    instanceSize = 0
    extraSize = 0
    # print(match)
    if match != None:
        instanceSize = int(match[0][0])
        extraSize = int(match[0][1])
    # print(instanceSize) [('200', '8')]

    start = time.time()
    # Create a graph, each node is a one k-gram
    content = []
    with open(argv[1]) as file:
        content = file.read().split('\n')[:-1]
    graph = nx.DiGraph()
    graph.add_nodes_from(content)

    maxPossibleLength = instanceSize + len(content[0]) - 1
    if "-" in argv[1]:
        extraSize = -extraSize
    else:
        extraSize = 0
    optimalWordsCount = instanceSize + extraSize

    # print("MAXLEN:", maxPossibleLength, "\nOPT:", optimalWordsCount, "\nFILE:", argv[1])
    
    edges = 0
    # Graph creation loop
    for seq in content:
            content2 = content[:]
            content2.remove(seq)
            for seq2 in content2:
                if ((overlap_count := check_overlap(seq, seq2)) > 1 ):
                    graph.add_edge(seq, seq2, weight = overlap_count)
                    edges += 1


    # Pruning the graph
    removedEdges = trimGraph(graph, aboveTreshPercentage=0.8, percentileLower=50, percentileUpper=90)
    
    # print(edges, removedEdges, edges-removedEdges)

    # Jak chodzimy:
    # 1. Po najlepszych wagach
    # 2. Po najlepszych następnikach
    #    - Najlepsza średnia wag u następnika
    #    - Najlepsza waga u następnika
    #    - Najwięcej następników u następnika
    #    - Najwięcej następników których jeszcze nie było
    # Jak długo chodzimy:
    # 1. Aż nie ma pętli
    # 2. Aż nie przekroczymy maksymalnej długości
    # 3. Aż nie osiągniemy OPT unikatowych wierzchołków

    # Plan:
    # Dla każdego wierzchołka:
    #     dodajemy wierzchołek do ścieżki
    #     idziemy do następnego wierchołka (2.4)
    #     wracamy do: dodajemy wierzchołek do ścieżki aż (2., 3.)
    #     zapisujemy wynik
    # Porównujemy wyniki bierzemy najbliższy od OPT
    # Print: Czas, odległość od OPT, Sekwencja

    bestPath, bestLen, distFromOpt = heuristicSucc(graph, optimalWordsCount, maxPossibleLength, 0)
    print(
"--------\n\
Filename: {file}\n\
OptLen = {oL} BestLen = {bL} DistFromOpt = {dFO}\n\
Best = {b}\n\
BestPathLen = {bPL}/{mPL}".format(
          file=argv[1],
          oL=optimalWordsCount,
          bL=bestLen,
          dFO=distFromOpt,
          b=bestPath,
          bPL=len(concatenatePath(bestPath)),
          mPL=maxPossibleLength))

    stop = time.time()
    print("Execution time:", str((stop - start)*1000) + "ms\n--------\n")

    graph.clear()