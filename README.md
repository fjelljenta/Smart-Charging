# Smart-Charging
Repository for our contribution to the Qiskit Global Hackathon 2021

During the hackathon we want to look at quantum algorithms for industrial NP-hard problems. Due to climate change we will focus on the impactful problem of smart charging. 


## The issue
Climate change is one of the greatest challenges of our time. In order to meet the Paris Agreement, all countries have to cut down their emissions drastically. Electric mobility will play a key role in order to reduce greenhouse gas emission in the transport sector. If more and more electric vehicles enter the streets, it will become more challenging to manage the power supply. Smart charging will be necessary to improve the flexibility of the electric grid. However, smart charging faces large sized combinatorial optimization problems and many of them are NP-hard. 

For our project we mainly followed a paper dealing with a case study in the field of smart charging of electrical vehicles by Constantin Dalyac and Loïc Henriet et al. (https://arxiv.org/pdf/2012.14859.pdf)

## What is Smart Charging?
Smart charging stands for the intelligent charging of electric cars and it has many facets. Electric vehicles should, for example, be charged when it makes economic sense. This includes charging operations at times when electricity is cheap and/or during load peaks, which means there is enough capacity available in the power grid. One also can take into account the needs of the car owner and use predictions when the car should be fully charged. Additionally, load management ensures that electric vehicles are optimally charged at their current parking location. Whether at home or at work – intelligent charging avoids high power grid upgrading costs and peak energy prices. Therefore, smart charging uses the resources optimally and makes it more comfortably for the car owner.

In the first part of this project, we focus on a timing problem. This means we want to optimize the time a set of cars needs to be charged. We imagine a charging station with a various number of charging points (k) and more cars/loading jobs (n) than available charging points. We now want to find the optimal order to charge these cars. To make the task more realistic, we create vehicles with high priority status, like ambulances or police cars, and cars with lower priority (w). Moreover, the vehicles have different charging times (t). In the end, we want to optimize the total weighted charging completion time for all vehicles in the set.

In the second part, we focus on an interval scheduling problem. Imaging a charging point which is shared by x companies. Each company expresses wishes for certain time intervals (start time, end time) during which they would like to have access to the loading station. We consider the case where all companies ask for the same number of time intervalls. The goal is to include as many companies as possible, so that no company is overrepresented and at the same time the completion time is minimized.  

## Minimization of Total Weighted Load Completion Time 
We consider  n electric vehicles, so we have J = {1, . . . , n} charging jobs with durations T = {t1, . . . , tn}. These charging jobs need to be scheduled on a set I = {1, . . . , k}, where k is the available number of charging points. To each job we assign a weight wj > 0, which denotes its importance. The time at which a load j ends, is called completion time Cj. Our goal is to minimize the weighted total time of completion (sum over Cj * wj), which we also call the cost of the problem.
This is a NP-hard scheduling problem. One can see this problem as a weighted Max-k-Cut problem, which can for example be tackled by a Quantum Approximation Optimization Algorithm (QAOA). 
For this purpose, we will translate the cost appearing in the Max-k-Cut problem into a cost Hamiltonian. For this we can find the ground state energy using a combined loop of quantum computing combined and classical optimization.

## Optimal Scheduling of Load Time Intervals within Groups
An optimal solution for the interval scheduling problem can be interpreted as a Maximal Independent Set. A Maximum Independent Set (MIS) of a Graph G is an indepent subset of nodes, such that no two nodes, appearing in the subset, are directly connected by an edge. Each node of the graph represents a loading job and has three attributes: the group (company) it belongs to, a start time and an end time of the charging process. There is a connection between two nodes if the two nodes are in the same group or/and the charging intervals are overlapping. The goal is to find a maximal subset of not connected nodes. Following the example in the paper, we will restric ourselves to MIS problems, which can be formulated on two dimensional unit disk graphes. In such graphes only these nodes are connected, which are close to each other quantified by a certain distant parameter. Again we can translate the search for the MIS into calculating the minimum of a cost function. Similary to the description for the Max-k-Cut problem, we can encode this cost function into a Hamiltonian and search for its lowest energy state.

## Overview over the repository
If you are now interested to have a deeper look into what we did in our project, we will recommend that you start by looking into the two Jupyter Notebooks. You can find one notebook for the Max-k-Cut problem including a benchmarking and one notebook for the MIS problem.
The functions we call in these notebooks are stored in the corresponding python files in the "Code" folder. In the "create_data" folder we stored data for most of the benchmark runs.
Furthermore, we stored our resources in the file "Resources.txt" if you wish to further read about the topics and implementations.

## Visualization
We additionally created a Django web app with a visual simulation of what our code does. Check it out if you are interested!
https://evchargingconsole.herokuapp.com/
You can find the code for the web app in in the folder "web_app/console_simulation".
