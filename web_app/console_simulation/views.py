from django.shortcuts import render
from django.http import JsonResponse
from . import mc_functions
from . import mis_functions
import json
from django.views.decorators.csrf import csrf_exempt


# Create your views here.

def index(request):
    return render(request, "console_simulation/index.html")

def maxcut_sim(request):
    return render(request, "console_simulation/maxcut_sim.html")

def mis_sim(request):
    return render(request, "console_simulation/mis_sim.html")

@csrf_exempt
def init_mc_jobs(request):
    post_data = json.loads(request.body.decode("utf-8"))
    n = int(post_data["n"])
    return JsonResponse(mc_functions.init_vehicles(n))

@csrf_exempt
def init_mis_jobs(request):
    post_data = json.loads(request.body.decode("utf-8"))
    n = list(post_data["n"])
    return JsonResponse(mis_functions.init_jobs(n))

@csrf_exempt
def get_maxcut(request):
    post_data = json.loads(request.body.decode("utf-8"))
    j = post_data["jobs"]
    # print(j)
    cut=mc_functions.get_mc(j)
    job_sequence={"p0":[],"p1":[],"p2":[],"p3":[]}
    for i in range(len(cut)):
        job_sequence["p"+str(cut[i])]+=[j[i]]
    for key in job_sequence.keys():
        print(job_sequence[key])
        job_sequence[key].sort(key=lambda job : float(job["charging_time"])/float(job["importance"]))
    return JsonResponse(job_sequence)

@csrf_exempt
def get_mis(request):
    jobs = json.loads(request.body.decode("utf-8"))
    return JsonResponse(mis_functions.get_mis(jobs))
