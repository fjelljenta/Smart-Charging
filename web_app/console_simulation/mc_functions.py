import random

def init_vehicles(n):
    d={"vehicles": [{"importance": round(0.1+0.9*random.random(), 3), "charging_time": round(5*(0.1+0.9*random.random()), 2)} for _ in range(n)]}
    return d

def get_mc(jobs):
    return [int(4*random.random()) for _ in jobs]
