"""
Assignment 7, Annotation-Based Pairwise Semantic Similarity Measure
2016253072 명수환
specificity of a term C_i : I(C_i)
"""


from datetime import datetime
from os import name
import sys
import math
import matplotlib.pyplot as plt
import threading
import multiprocessing
import time

def getDataFromFile(filename):
    file = open(filename, "r")

    ontology_MF = dict()
    ontology_BP = dict()

    line = None
    line = file.readline()
    # print(line)
    while line != "":
        # print(line)
        while line != "[Term]\n":
            line = file.readline()  # * Term 단위
            if line == "":
                return ontology_MF, ontology_BP

        line = file.readline().split()
        # print(line)
        # * Read line after "[Term]\n"
        if line[0] == "id:":
            id = line[1]
            # print("id:",id)
            line = file.readline().split()

        if line[0] == "name:":
            # print("name:",line[1:])
            line = file.readline().split()

        if line[0] == "namespace:":
            namespace = line[1]
            if namespace != "molecular_function" and namespace != "biological_process":
                line = file.readline().split()
                continue
            # print("namespace:",namespace)#* molecular_function,biological_process
            line = file.readline().split()

        is_a_relation = set()
        part_of_relation = set()
        is_obsleted = False
        while line:
            # print(line)
            if line[0] == "is_obsolete:":
                if line[1] == "true":
                    # print(id,"is obsolete")
                    line = file.readline().split()
                    is_obsleted = True
                    break
            elif line[0] == "is_a:":
                is_a_relation.add(line[1])
                line = file.readline().split()
            elif line[0] == "relationship:":
                if line[1] == "part_of":
                    part_of_relation.add(line[2])  # * id
                    line = file.readline().split()
                else:
                    line = file.readline().split()
            else:
                line = file.readline().split()

        if is_obsleted:
            continue

        error_relation = is_a_relation.intersection(part_of_relation)
        if error_relation:
            print(id, "has two relationship - ", error_relation)

        # * Classify - BP,MF
        if namespace == "molecular_function":
            ontology_MF[id] = list(
                is_a_relation.union(part_of_relation)
            )  # * {id:relation}
        elif namespace == "biological_process":
            ontology_BP[id] = list(is_a_relation.union(part_of_relation))

    return ontology_MF, ontology_BP


def getTermDepth(ontology, start_node, root_node):
    shortestPathLength = -1  # * 초기값

    for v in ontology[start_node]:
        edges_count = 0
        parent = v
        # print("parent : ", parent)
        edges_count += 1
        if parent != root_node:
            start_node = parent
            edges_count += getTermDepth(ontology, start_node, root_node)
        else:  # * parent is root
            return 1
        # print("short: ",shortestPathLength)
        # print("edge count : ",edges_count)
        if shortestPathLength != -1:
            if edges_count < shortestPathLength:
                shortestPathLength = edges_count
        elif shortestPathLength == -1:
            shortestPathLength = edges_count
        else:  # * shortestPathLength < edges_count
            pass

    # print("short: ",shortestPathLength)
    return shortestPathLength


def getDataFromGAF(filename, BP, MF):

    file = open(filename, "r")

    line = None

    while True:

        line = file.readline()

        if not line:
            break

        if line[0] == "!":
            continue

        line = line.split("\t")

        term = line[4].strip()
        gene = line[2].strip()
        code = line[6].strip()
        confirm = line[3].strip()
        ontology = line[8].strip()

        if gene == "":

            continue
        if term == "":

            continue
        if code == "":

            continue
        if code == "IEA":

            continue
        if confirm[:3] == "NOT":

            continue
        if ontology == "C":

            continue

        if ontology == "F":

            MF[term].add(gene)

        elif ontology == "P":
            BP[term].add(gene)

    return BP, MF


def GetPPI(file, annotation_BP, annotation_MF, ontology_BP, ontology_MF,start_position):
    global global_gene_sim
    
    #fileopen
    line_count = 0
    total_line = 159567
    line = None
    gene_sim = list()

    # * for starting position
    for __ in range(start_position):
            file.readline()
    print("Thread ", start_position, "is started at position ",file.tell())
    while True:
        
        line_count += 1
        if line_count % int(total_line / 100) == 0:
            print("Progress : {0}%".format(round(line_count * 100 / total_line, 3)))
            # print("line : ", line)
        line = file.readline().split()
        
        # * for threading
        for _ in range(10-1):
            file.readline()

        print("line : ", line)
        if not line:
            break
        gene1 = line[0]
        gene2 = line[1]

        gene1_BP = dict()
        gene2_BP = dict()
        gene1_MF = dict()
        gene2_MF = dict()
        # * find all terms having gene1 or gene2
        for term in annotation_BP:
            if gene1 in annotation_BP[term]:
                if not gene1 in gene1_BP:
                    gene1_BP[gene1] = set()

                gene1_BP[gene1].add(term)

            if gene2 in annotation_BP[term]:
                if not gene2 in gene2_BP:
                    gene2_BP[gene2] = set()
                gene2_BP[gene2].add(term)

        for term in annotation_MF:

            if gene1 in annotation_MF[term]:
                if not gene1 in gene1_MF:
                    gene1_MF[gene1] = set()
                gene1_MF[gene1].add(term)

            if gene2 in annotation_MF[term]:
                if not gene2 in gene2_MF:
                    gene2_MF[gene2] = set()
                gene2_MF[gene2].add(term)
        # * inferred 방식이면 위에서 이미 해당 gene을 갖는 모든 노드가 구해짐.
        # * GetAllAncestor 함수 호출할 필요가 X
        # * gene1 -> 구해진 모든 노드들 중에서 maximum I(C)를 갖는 노드를 구해서 C1
        # * gene2 -> C2
        # print("gene1 : ", gene1)
        # print("gene2 : ", gene2)
        C1_BP = set()
        C2_BP = set()
        C1_MF = set()
        C2_MF = set()

        sim_BP_list = list()
        sim_MF_list = list()
        max_MF_list = list()
        max_BP_list = list()

        if gene1 in gene1_BP and gene2 in gene2_BP:
            # print("BP case ")
            for v1bp in gene1_BP[gene1]:
                C1_BP = v1bp
                for v2bp in gene2_BP[gene2]:

                    C2_BP = v2bp
                    sim_BP = GetSimilarity_annotation_based(
                        C1_BP, C2_BP, ontology_BP, annotation_BP
                    )

                    # sim_BP_list.append([C1_BP, C2_BP, sim_BP])
                    sim_BP_list.append(sim_BP)
                    # print(sim_BP_list)
                max_BP_list.append(max(sim_BP_list))
            X_len_BP = len(max_BP_list)

            for v1bp in gene2_BP[gene2]:
                C1_BP = v1bp
                for v2bp in gene1_BP[gene1]:

                    C2_BP = v2bp
                    sim_BP = GetSimilarity_annotation_based(
                        C1_BP, C2_BP, ontology_BP, annotation_BP
                    )

                    # sim_BP_list.append([C1_BP, C2_BP, sim_BP])
                    sim_BP_list.append(sim_BP)
                max_BP_list.append(max(sim_BP_list))
            Y_len_BP = len(max_BP_list) - X_len_BP

        if gene1 in gene1_MF and gene2 in gene2_MF:
            # print("MF case ")
            # * gene1을 갖는 term하나, gene2를 갖는 term하나를 구해서
            # * annotation_based method로 similarity 계산
            for v1mf in gene1_MF[gene1]:
                C1_MF = v1mf
                for v2mf in gene2_MF[gene2]:

                    C2_MF = v2mf
                    sim_MF = GetSimilarity_annotation_based(
                        C1_MF, C2_MF, ontology_MF, annotation_MF
                    )

                    sim_MF_list.append(sim_MF)
                max_MF_list.append(max(sim_MF_list))
            X_len_MF = len(max_MF_list)
            # * 여기서 max값 구하기

            for v1mf in gene2_MF[gene2]:
                C1_MF = v1mf
                for v2mf in gene1_MF[gene1]:

                    C2_MF = v2mf
                    sim_MF = GetSimilarity_annotation_based(
                        C1_MF, C2_MF, ontology_MF, annotation_MF
                    )

                    sim_MF_list.append(sim_MF)
                max_MF_list.append(max(sim_MF_list))
            Y_len_MF = len(max_MF_list) - X_len_MF

        if sim_BP_list and sim_MF_list:
            # * 둘다
            # max_similarity_BP = max(sim_BP_list)
            # max_similarity_MF = max(sim_MF_list)
            # max_similarity = max(max_similarity_BP,max_similarity_MF)

            BP_BMA = BestMatchingAvg(max_BP_list, X_len_BP, Y_len_BP)
            MF_BMA = BestMatchingAvg(max_MF_list, X_len_MF, Y_len_MF)
            if BP_BMA > MF_BMA:
                BMA = BP_BMA
            else:
                BMA = MF_BMA

        elif sim_BP_list and not sim_MF_list:
            BP_BMA = BestMatchingAvg(max_BP_list, X_len_BP, Y_len_BP)
            BMA = BP_BMA

        elif sim_MF_list and not sim_BP_list:
            MF_BMA = BestMatchingAvg(max_MF_list, X_len_MF, Y_len_MF)
            BMA = MF_BMA
        else:
            # print("error case")
            continue

        similarity = BMA
        gene_sim.append([gene1, gene2, similarity])
        # * gene_sim += node_based_method(gene1,gene2,gene1_BP,gene2_BP,gene1_MF,gene2_MF,ontology_BP,ontology_MF)

    global_gene_sim += gene_sim
    #return gene_sim


def BestMatchingAvg(max_similarity_list, X_len, Y_len):
    # print(" len of list ", len(max_similarity_list))
    # print(" Xlen+Ylen : ", (X_len+Y_len))
    similarity = sum(max_similarity_list) / (X_len + Y_len)
    # print("max_sim : ",max_similarity)
    # *print("BMA similarity : ", similarity)
    return similarity


def node_based_method(
    gene1, gene2, gene1_BP, gene2_BP, gene1_MF, gene2_MF, ontology_BP, ontology_MF
):
    # * Node based method, groupwise
    gene_sim = list()

    gene1_BP_ancestor = set()
    gene2_BP_ancestor = set()
    gene1_MF_ancestor = set()
    gene2_MF_ancestor = set()

    if gene1_BP:

        gene1_BP_ancestor = GetAllAncestor(gene1_BP[gene1], ontology_BP)

    if gene2_BP:
        gene2_BP_ancestor = GetAllAncestor(gene2_BP[gene2], ontology_BP)
    if gene1_MF:
        gene1_MF_ancestor = GetAllAncestor(gene1_MF[gene1], ontology_MF)
    if gene2_MF:
        gene2_MF_ancestor = GetAllAncestor(gene2_MF[gene2], ontology_MF)

    if (
        gene1_BP_ancestor
        and gene1_MF_ancestor
        and gene2_BP_ancestor
        and gene2_MF_ancestor
    ):
        # * Larger case within BP and MF
        # print("gene MF or BP")
        # print(gene1_BP_ancestor)
        # print(gene2_BP_ancestor)
        # print(gene1_MF_ancestor)
        # print(gene2_MF_ancestor)
        sim_BP = GetSimilarity_node_based(gene1_BP_ancestor, gene2_BP_ancestor)
        sim_MF = GetSimilarity_node_based(gene1_MF_ancestor, gene2_MF_ancestor)
        sim = max(sim_BP, sim_MF)
        tmp_list = [gene1, gene2, sim]
        gene_sim.append(tmp_list)
        # print("sim : ",sim)

    elif gene1_MF_ancestor and gene2_MF_ancestor:
        # * MF
        # print("gene MF")
        # print("gene1_MF_ancestor : ",gene1_MF_ancestor)
        # print("gene2_MF_ancestor : ",gene2_MF_ancestor)
        sim_MF = GetSimilarity_node_based(gene1_MF_ancestor, gene2_MF_ancestor)
        tmp_list = [gene1, gene2, sim_MF]
        gene_sim.append(tmp_list)
        # print("sim : ",sim_MF)

    elif gene1_BP_ancestor and gene2_BP_ancestor:
        # * BP
        sim_BP = GetSimilarity_node_based(gene1_BP_ancestor, gene2_BP_ancestor)
        tmp_list = [gene1, gene2, sim_BP]
        gene_sim.append(tmp_list)
        # print("sim : ",sim_BP)
    else:
        gene_sim = list()
        return gene_sim

        # * BP and MF
        # * or MF and BP
        # print("NOT Case!")
    # print("gene_sim = ", gene_sim)
    return gene_sim


def GetAllAncestor(nodeset, ontology):
    if not nodeset:
        return set()
    newset = set()

    for node in nodeset:

        if ontology[node]:
            for v in ontology[node]:
                newset.add(v)

    newset = newset.union(GetAllAncestor(newset, ontology))

    return newset


def GetCommonAncestor(term1, term2, ontology):
    # print("term1 : ", term1)
    # print("term2 : ", term2)
    set1 = {term1}
    set2 = {term2}
    ancestor1 = GetAllAncestor(set1, ontology)
    ancestor2 = GetAllAncestor(set2, ontology)
    # print("ancestor1 : ", ancestor1)
    # print("ancestor2 : ", ancestor2)
    if not ancestor1:
        return set()
    if not ancestor2:
        return set()
    common = ancestor1.intersection(ancestor2)

    #! term중 하나가 루트인 경우, common ancestor가 없으면 어떡하지?
    return common


def GetSimilarity_node_based(gene1_ancestor, gene2_ancestor):
    union_set = gene1_ancestor.union(gene2_ancestor)

    inter_set = gene1_ancestor.intersection(gene2_ancestor)
    # print("union : ",union_set)
    # print("inter : ",inter_set)
    return round(float(len(inter_set)) / float(len(union_set)), 1)


def GetSimilarity_annotation_based(C1, C2, ontology, annotation):
    # * sim(c1,c2) = 2*math.log(P(C0)) / (math.log(P(C1)) + math.log(P(C2)))
    # print("GetSImilarity func ..  .")
    P_C1 = P(C1, annotation)
    P_C2 = P(C2, annotation)
    I_term = 0
    common_ancestors = GetCommonAncestor(C1, C2, ontology)
    if not common_ancestors:
        return 0
    for ancestor in common_ancestors:
        # print("ancestor : ", ancestor)
        tmp = GetInformationContent(ancestor, annotation)
        # print(tmp)
        if I_term <= tmp:  # * if I_term is -0.0 (root case)
            I_term = tmp
            # print("I_term : ", I_term)
            C0 = ancestor
            # print("C0 : ", ancestor)
    P_C0 = P(C0, annotation)
    # * P_C0 : maximum information content. (the information content of the most specific common ancestor
    similarity = 2 * math.log(P_C0) / (math.log(P_C1) + math.log(P_C2))
    # print("similarity : ", similarity)
    # * information concept : P(C) = C에 anootated된 gene의 수  / 온톨로지 전체에 annotated 된 genes의 수
    return similarity


def output_to_file(filename, gene_sim):
    file = open(filename, "w")

    for v in gene_sim:
        file.write("{0} {1} : {2}".format(v[0], v[1], v[2]))
        file.write("\n")
    file.close()


def GetInformationContent(term, anotation_):
    # * I(C) = -logP(C)

    P_term = P(term, anotation_)
    # print("P_term : ", P_term)
    I_term = (-1) * math.log(P_term)
    # print("I term : ",I_term)

    return I_term


def P(term, annotation_):
    if root_BP in annotation_:
        total_gene_num = BP_gene_set_length
    elif root_MF in annotation_:
        total_gene_num = MF_gene_set_length
    
    P_term = len(annotation_[term]) / total_gene_num

    return P_term

def inferred(ontology,annotation,total_length,root_node):
    print("[*]inferred start")
    while True:
        if len(annotation[root_node]) == total_length:
            satisfied = True
        else:
            satisfied = False
        

        if satisfied:
            break

        if not satisfied:
            for child in annotation:
                for parent in ontology[child]:
                    annotation[parent] = annotation[parent].union(
                        annotation[child]
                    )
    
    print("[*]inferred end")    


def GetTotalLengthOfGenes(annotation,root_node,root_BP,root_MF):
    

    gene_set_ = set()
    for term in annotation:
        gene_set_ = gene_set_.union(annotation[term])
    
    if root_node == root_BP:
        global BP_gene_set_length
        BP_gene_set_length = len(gene_set_)
    if root_node == root_MF:
        global MF_gene_set_length
        MF_gene_set_length = len(gene_set_)

def main():

    if len(sys.argv) != 4:
        print("No input file.")
        print("<Usage> Assignment7.py go.obo goa_human.gaf biogrid_human_ppi_cln.txt")
        return -1

    input_filename = sys.argv[1]
    input_filename2 = sys.argv[2]
    input_filename3 = sys.argv[3]
    output_filename = "output.txt"

    # input_filename = "go.obo"
    # input_filename2 = "goa_human.gaf"
    # input_filename3 = "biogrid_human_ppi_cln.txt"

    start_time = datetime.now()

    ontology_MF, ontology_BP = getDataFromFile(input_filename)

    print("length of MF :", len(ontology_MF))
    print("length of BP :", len(ontology_BP))

    global root_MF
    global root_BP

    for v1 in ontology_MF:

        if ontology_MF[v1] == []:
            root_MF = v1
            print("root node in MF : ", v1)

    for v2 in ontology_BP:

        if ontology_BP[v2] == []:
            root_BP = v2
            print("root node in BP : ", v2)

    error_count = 0
    for id_bp in ontology_BP:

        if id_bp == root_BP:
            continue
        for u1 in ontology_BP[id_bp][:]:

            if u1 in ontology_MF:
                # * ERROR -  Relation u and  v :  u is in MF, v is in BP
                error_count += 1
                ontology_BP[id_bp].remove(u1)
                # print(id_bp,"에서 ",u1,"를 삭제하였습니다.")

    for id_mf in ontology_MF:

        if id_mf == root_MF:
            continue
        for u2 in ontology_MF[id_mf][:]:
            if u2 in ontology_BP:
                # * ERROR - Relation u and  v  :  u is in BP, v is in MF
                error_count += 1
                ontology_MF[id_mf].remove(u2)

                # print(id_mf,"에서 ",u2,"를 삭제하였습니다.")
    print("error count : ", error_count)

    BP_annotation = ontology_BP.copy()
    MF_annotation = ontology_MF.copy()
    for v in BP_annotation:
        BP_annotation[v] = set()
    for v in MF_annotation:
        MF_annotation[v] = set()  # * label이 없는 term이 존재할 수 있음.

    BP_annotation, MF_annotation = getDataFromGAF(
        input_filename2, BP_annotation, MF_annotation
    )

    
    print("=========================================================================")
    print("length of BP_annotation : ", len(BP_annotation))
    print("length of MF_annotation : ", len(MF_annotation))
    global MF_gene_set_length
    global BP_gene_set_length
    
    
    threads = []
    t1 = threading.Thread(target=GetTotalLengthOfGenes,args=(MF_annotation,root_MF,root_BP,root_MF))
    t1.start()
    threads.append(t1)
    
    t2 = threading.Thread(target=GetTotalLengthOfGenes,args=(BP_annotation,root_BP,root_BP,root_MF))
    t2.start()
    threads.append(t2)

    for thread in threads:
        thread.join()
    


    # process_ = []
    # p1 = multiprocessing.Process(target=GetTotalLengthOfGenes,args=(MF_annotation,root_MF,root_BP,root_MF))
    # p2 = multiprocessing.Process(target=GetTotalLengthOfGenes,args=(BP_annotation,root_BP,root_BP,root_MF))

    # p1.start()
    # p2.start()

    # process_.append(p1)
    # process_.append(p2)

    # for proc in process_:
    #     proc.join()
    
    print("MF_gene_set_length : ", MF_gene_set_length)
    print("BP_gene_set_length : ", BP_gene_set_length)

    print("BEFORE root_MF labels length ", len(MF_annotation[root_MF]))
    print("BEFORE root_BP labels length ", len(BP_annotation[root_BP]))


    # * inferred
    threads = []
    t1 = threading.Thread(target=inferred,args=(ontology_BP,BP_annotation,BP_gene_set_length,root_BP))
    t1.start()
    threads.append(t1)
    
    t2 = threading.Thread(target=inferred,args=(ontology_MF,MF_annotation,MF_gene_set_length,root_MF))
    t2.start()
    threads.append(t2)

    for thread in threads:
        thread.join()
    print("end inferred")

    # * show annotation

    print("AFTER root_MF labels length ", len(MF_annotation[root_MF]))
    print("AFTER root_BP labels length ", len(BP_annotation[root_BP]))


    global global_gene_sim
    global_gene_sim = list()

    
    file_PPI = open(input_filename3, "r")
    threads = []
    print("PPI READ start")
    for start_pos in range(10):
        print("Thread ", start_pos, " START ")
        t = threading.Thread(target=GetPPI,args=(file_PPI, BP_annotation, MF_annotation, ontology_BP, ontology_MF, start_pos))
        t.start()
        threads.append(t)
    
    for thread in threads:
        thread.join()
    print("PPI READ end")
    
    
    # process_ = []
    # for start_pos in range(10):
    #     print("Process ", start_pos, " START ")
    #     p = multiprocessing.Process(target=GetPPI, args=(file_PPI, BP_annotation, MF_annotation, ontology_BP, ontology_MF, start_pos))
    #     p.start()
    #     process_.append(p)
    
    # for proc in process_:
    #     proc.join()
    
    file_PPI.close()

    similarity = list()
    for v in global_gene_sim:
        similarity.append(v[2])
    print("Similarity 0.0 ~ 1.0 : ", similarity)
    print("[+] Time Elapsed : ", datetime.now() - start_time, "microseconds")
    plt.hist(similarity)
    plt.show()
    # output_to_file(output_filename,gene_sim)


if __name__ == "__main__":
    main()
