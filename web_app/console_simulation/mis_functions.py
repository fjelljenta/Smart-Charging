import random
from console_simulation.Maximal_independent_set_quantum_functions import *

def init_jobs(n):
    d={}
    for company in range(len(n)):
        d["company-"+chr(ord('a') + company)]=[{"start_time": 7+8*random.random(), "duration": 1+3*random.random()} for job in range(n[company])]
    return d

def get_mis(jobs):
    # mis=[]
    # for company in jobs.keys():
    #     if random.random()>0.3:
    #         index=1+int(len(jobs[company])*random.random())
    #         mis+=[{"job-ID": "job-"+company[-1].upper()+str(index), "start_time": jobs[company][index-1]["start_time"]}]
    # mis.sort(key = lambda el : el["start_time"])

    p = 2  # the depth of QAOA
    ng = 4  # number_of_groups
    n = 10  # total number of nodes, or vehicles
    mt = 500  # max_time, Start time charging limit (say, after 100 minutes, every car started at least charging)
    mct = 100  # max_charging_time, maximum time a car can charge (like 6h)
    U = 10  # The size of punishment added to the cost due to including nodes of the same group in the MIS

    j=[]

    for ch in ['a','b','c','d']:
        j+=[{**job, **{"company": "company-"+ch, "index": index}} for index, job in enumerate(jobs["company-"+ch])]

    adj = []

    for i in range(10):
        row=[]
        for k in range(10):
            row += [1] if i!=k and (j[i]["company"]==j[k]["company"] or min(float(j[i]["start_time"])+float(j[i]["duration"]), float(j[k]["start_time"])+float(j[k]["duration"]))>max(float(j[i]["start_time"]), float(j[k]["start_time"]))) else [0]
        adj+=[row]

    print(adj)

    circ = make_full_circuit_MIS(n,adj,p)
    simulator='qasm_simulator'
    backend = Aer.get_backend(simulator)
    backend.shots = 128
    circ_w_param = circ.bind_parameters([ 0.77465785, 2.82478591, -1.13516493, 3.90552469])
    transpiled_circ = transpile(circ_w_param, backend)
    counts = backend.run(transpiled_circ, shots=128).result().get_counts()
    measurement=max(counts, key=counts.get)

    preprocessed_chosen_set = np.argwhere(np.array(list(measurement))=='1')
    chosen_set = np.concatenate(preprocessed_chosen_set)

    return {"mis": ["job-"+j[pos]["company"][-1].upper()+str(j[pos]["index"]+1) for pos in chosen_set]}

    # return {"mis": [el["job-ID"] for el in mis]}
