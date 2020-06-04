from sys import argv
import networkx as nx
import matplotlib.pyplot as plt

def check_overlap(seq: str, seq2: str):
    init_overlap = len(seq)
    for i in range(len(seq2)):
        if seq.endswith(seq2[:init_overlap]):
            return init_overlap
        init_overlap -= 1
    
    return init_overlap
    




if __name__ == "__main__":
    if (len(argv) < 2):
        print("Usage: " + argv[0] + " [filename]")
        exit(1)
    
    content = []
    with open(argv[1]) as file:
        content = file.read().split('\n')[:-1]
    graph = nx.DiGraph()
    graph.add_nodes_from(content)
    
    print(graph.number_of_nodes())
    

    for seq in content:
        content2 = content[:]
        content2.remove(seq)
        for seq2 in content2:
            if ((overlap_count := check_overlap(seq, seq2)) > 6):
                graph.add_edge(seq, seq2)
                print(seq, seq2)
                graph.add_weighted_edges_from([(seq,seq2,overlap_count)])
                #graph[seq][seq2]['weight'] = overlap_count
                #print(seq, seq2, overlap_count)
            #TODO: check overlap and build a overlap graph

    for e in graph.edges(data='weight'):
        print(e)
        #print(e, graph[e[0]][e[1]]['weight'])
        

    nx.draw(graph)
    plt.show()

    print(content)
    graph.clear()