import random

def init_jobs(n):
    d={}
    for company in range(len(n)):
        d["company-"+chr(ord('a') + company)]=[{"start_time": 7+8*random.random(), "duration": 1+3*random.random()} for job in range(n[company])]
    return d

def get_mis(jobs):
    mis=[]
    for company in jobs.keys():
        if random.random()>0.3:
            index=1+int(len(jobs[company])*random.random())
            mis+=[{"job-ID": "job-"+company[-1].upper()+str(index), "start_time": jobs[company][index-1]["start_time"]}]
    mis.sort(key = lambda el : el["start_time"])
    return {"mis": [el["job-ID"] for el in mis]}
