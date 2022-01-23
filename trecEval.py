# We implement in this file all the evaluation metrics used in our experiments including MAP, P@k, nDCG and nDCG@k.

import numpy as np
import pandas as pd
import sys


def log2(x):
    from math import log
    return log(x, 2)


def getQrels(qrelsFile):
    with open(qrelsFile, 'r') as f:
        qrels = {}
        QIDs = []
        for line in f:
            content = line.split()
            qID = content[0]
            docID = content[2]
            rel = content[3]
            if qID not in qrels:
                QIDs.append(qID) 
                qrels[qID] = {}
                docs = qrels[qID]
            if int(rel):
                docs[docID] = rel
    return QIDs, qrels


def getResults(resultsFile):
    with open(resultsFile, 'r') as f:
        results = {}
        for line in f:
            content = line.split()
            qID = content[0]
            docID = content[2]
            if qID not in results:
                results[qID] = [docID]
            else:
                temp = results[qID]
                temp.append(docID)
                results[qID] = temp
    return results


def AP(rels, ranked_list):
    ap = 0
    found = 0
    if not rels:
        return np.nan
    for i, docID in enumerate(ranked_list):
        if docID in rels:
            found += 1
            ap += found / (i+1)
    ap = ap/len(rels)
    return ap


def prec(rels, ranked_list, k):
    p = 0
    found = 0
    if not rels:
        return np.nan
    for i, docID in enumerate(ranked_list):
        if i < k:
            if docID in rels:
                found += 1
    p = found/k
    return p


def NDCG(rels, ranked_list, k):
    dcg = 0
    ndcg = 0
    n = 0
    found = 0
    if not rels:
        return np.nan
    for i, docID in enumerate(ranked_list):
        if (i < k) and (i < len(rels)) and (found < len(rels)):
            if docID in rels:
                found += 1
                dcg += 1.0/log2(i+2.0)
    for i, docID in enumerate(rels):
        if i < k:
            n += 1.0/log2(i+2.0)
    ndcg = dcg/n
    return ndcg


# res is a map{qID, value}
def Mean(res):
    avg = 0
    i = 0
    for qID, value in res.items():
        if not np.isnan(value):
            i += 1
            avg += value
    avg = avg/i
    return avg


# res is list [value]
def average(res):
    avg = 0
    i = 0
    for value in res:
        if not np.isnan(float(value)):
            i += 1
            avg += float(value)
    avg = avg/i
    return avg


def eval_one_file(resultFile, qrelsFile, outFile):
    QIDs, qrels = getQrels(qrelsFile)

    results = getResults(resultFile)
    temp = list(qrels.items())

    for QID, value in temp:
        if not value:
            QIDs.remove(QID)
            del qrels[QID]

    data = []
    for qID in QIDs:
        ap = round(AP(qrels[qID], results[qID]), 4)
        p5 = round(prec(qrels[qID], results[qID], 5), 4)
        p10 = round(prec(qrels[qID], results[qID], 10), 4)
        p15 = round(prec(qrels[qID], results[qID], 15), 4)
        p20 = round(prec(qrels[qID], results[qID], 20), 4)
        ndcg = round(NDCG(qrels[qID], results[qID], 1000), 4)
        ndcg5 = round(NDCG(qrels[qID], results[qID], 5), 4)
        ndcg10 = round(NDCG(qrels[qID], results[qID], 10), 4)
        ndcg15 = round(NDCG(qrels[qID], results[qID], 15), 4)
        ndcg20 = round(NDCG(qrels[qID], results[qID], 20), 4)
        temp = [qID, ap, p5, p10, p15, p20, ndcg, ndcg5, ndcg10, ndcg15, ndcg20]
        data.append(temp)
    temp = ["Average"]
    for i in range(1, len(data[0])):
        temp.append(round(average(np.array(data)[:, i]), 4))

    data.append(temp)
    df = pd.DataFrame(data,
                      columns=["qID", "AP", "P@5", "P@10", "P@15", "P@20",
                               "NDCG", "NDCG@5", "NDCG@10", "NDCG@15", "NDCG@20"],
                      index=None)
    writer = pd.ExcelWriter(outFile)
    df.to_excel(writer, index=False)
    writer.save()


if __name__ == '__main__':
    argv = sys.argv
    # "result/BM25/WT2G-BM25-1.2-0.35-report.txt"
    resultFile = argv[1]
    # "query-judge/qrels.WT2G"
    qrelsFile = argv[2]
    # "result/WT2G-BM25-1.2-0.35.xls"
    outFile = argv[3]
    eval_one_file(resultFile, qrelsFile, outFile)
    
