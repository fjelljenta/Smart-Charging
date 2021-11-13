# Smart-Charging
Repository for our contribution to the Qiskit Global Hackathon 2021

During the hackathon we wanted to look at quantum algorithms for industrial NP-hard problems. Due to climate change we will focus on the impactful problem of smart charging. 

## The issue
Climate change is one of the greatest challenges of our time. In order to meet the Paris Agreement, all countries have to cut down their emissions drastically. Electric mobility will play a key role in order to reduce greenhouse gas emission in the transport sector. If more and more electric vehicles enter the streets, it will become more challenging for the power grid. Smart charging will be necessary to improve the flexibility of the electric grid. However, smart charging faces large sized combinatorial optimization problems, many of them are NP-hard. 

## What is Smart Charging?
<<<<<<< HEAD
Smart Charging means the intelligent charging of electric cars. Smart charging has a many faccets. Electric vehicles should, for example, be charge when it makes economic sense. This includes charging operations at times when electricity is cheap but also charging it when it is possible, meaning there is enough capacity available on the power grid. One also can take in account the needs of the car owner and use predictions when the car should be fully charged. Additionally, load management ensures that electric vehicles are optimally charged at their location. Whether at home or at work – intelligent charging avoids high power grid upgrading costs, peak energy prices. Therefore using the resources optimally and making it more comfortably for the car owner.

In this project we focus on a timing problem. This means we want to optimize the time a set of cars need to be charged. We imagine a charging station with a various number of charging points and more cars than available charging points. We now want to find the optimal order to charge these cars. To make the task more realistic, we take vehicles with high priority status, like ambulances or police cars, with into account and give them a higher priority. In the end, we want the optimize the total weighted charging completion time for all vehicles in the set.
# In this project we focus on a problem regarding the car owner's perspective. We imagine a loading station, where more cars are than charging points. In this case one has to find an optimal way to charge the cars. We assume to have two attributes: the importance of a car (e.g if it is ambulance or police) and the forecast charging time. This results in the problem of the minimization of the total weighted load completion time.
=======
Smart Charging means the intelligent charging of electric cars. In the case of smart charging, electric vehicles always charge when it makes economic sense. On the one hand, this includes charging operations at times when electricity is cheap. On the other hand, smart charging also means load management: the adaptation of charging speed and performance to the available grid capacity and the needs of the car owner. Load management ensures that electric vehicles are optimally charged at their location. Whether at home or at work – intelligent charging avoids high network upgrading costs and peak prices and makes optimal use of the existing network connection.

In this project we focus on a problem regarding the car owner's needs. We imagine a situation, where we have more cars than charging points. In this case one has to find an optimal way to charge the cars. We assume to have two attributes: the importance of a car (e.g if it is ambulance or police) and the forecast charging time. This results in the problem of the minimization of the total weighted load completion time.
>>>>>>> 07e06740968c2522536b5754603d54182fc1b4c2

## Minimization of Total Weighted Load Completion Time 
We consider  n electric vehicles, so we have J = {1, . . . , n} charging jobs with durations T = {t1, . . . , tn}. These charging jobs need to be scheduled on a set I = {1, . . . , k}, where k is the number of charging points. To each job we assign a weight wj > 0, which gives the importance. The time at which a load j ends, is called completion time Cj. Our goal is to minimize the weighted total time of completion.
This is a classical NP-hard scheduling problem. One can see this problem as a weighted Max-k-Cut problem, which can be tackled by the Quantum Approximation Optimization Algorithm (QAOA), at least in its k = 2 set up i.e. Max-Cut.
