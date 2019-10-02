import pickle
import pandas


query = pickle.load(open("query1.pkl", "rb"))
result1 = pickle.load(open("result1_GEO.pkl", "rb"))
result2 = pickle.load(open("result1_MC.pkl", "rb"))


idx = [q[0] for q in query]
n = [len(q[-1]) for q in query]

import numpy as np

r1 = [result1[i][0] if result1.__contains__(i) else np.nan for i in idx]
r2 = [result2[i][0] if result2.__contains__(i) else np.nan for i in idx]


t1 = [result1[i][1] if result1.__contains__(i) else -1 for i in idx]
t2 = [result2[i][1] if result2.__contains__(i) else -1 for i in idx]

df = pandas.DataFrame({"index": idx, "number of fov": n, "result_GEO": r1, "time1": t1, "result_MC": r2, "time2": t2, "error": (np.array(r1) - np.array(r2))})
df.to_csv("exp_MC.csv")