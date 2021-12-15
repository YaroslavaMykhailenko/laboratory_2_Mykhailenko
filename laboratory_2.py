import random
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from seaborn.palettes import color_palette
from statistics import mean, variance
import mpmath as mp
import copy



#  1е завдання
def PuassonProcess():

  arr_sk = []    
  arr_time = []  
  time = 0

  lambda_ = float(input("Enter intensity: "))
  n = int(input("Enter number of events: "))

  for i in range(n):
    S_k = -math.log(1.0 - random.random()) / lambda_
    arr_sk.append(S_k)
    time = time + S_k
    arr_time.append(time)

#  графік реалізацій процесу

  plt.plot(np.arange(0, 1, 1.0 / n), arr_time)
  plt.ylabel('time')
  plt.show()


#  1.1 

  event_num = int(input("Enter number of the event: "))
  
  sk_tmp = arr_sk.copy()
  time_appear = []

  for i in range(len(sk_tmp)):
      element = sum(sk_tmp[:event_num])
      del sk_tmp[:event_num]
      time_appear.append(element)

#  гістограма розподілів часу появи заданої події 

  sns.histplot(time_appear,  edgecolor = "black", color= "pink", bins=10)
  plt.show()


 
 # 1.2

  interval_events = []
  for i in range(1,len(arr_sk)):
    interval_events.append(abs(arr_sk[i]-arr_sk[i-1]))

#  гістограма розподілів інтервалу між подіями  

  sns.histplot(interval_events,  edgecolor = "black", color= "pink",bins=10)
  plt.show()

# 1.3

  n = 10
  T = n/lambda_
  occur_event = []
  events = []
  tmp_sk_arr = arr_sk.copy()

  for i in range(int(len(tmp_sk_arr) / n)):
   sum_ = sum(tmp_sk_arr[:n])
   if sum_ < T and sum(tmp_sk_arr[:n+1]) > T:
    occur_event.append(1)
    events.append(sum_)
    del tmp_sk_arr[:n]

   else:
    occur_event.append(0)
    del tmp_sk_arr[:n]


#     # гістограма розподілів появи рівно n-подій 
  freq = occur_event.count(1)/len(occur_event)
  freq = (math.exp(-lambda_ * T) * pow((lambda_*  T),n)) / math.factorial(n)


  sns.histplot(events,  edgecolor = "black", color= "pink", bins=10)
  plt.show()











# 2е завдання
def WienerProcess():
    M = 1000
    T = 1
    t = T / M
    eta = np.random.normal(0, 0.1, M)
    eta1 = np.random.normal(0, 0.1, M)
    eta2 = np.random.normal(0, 0.1, M)

    def wiener1(eta):
        first_ksi_arr = []
        currenttime = 0

        for _ in range(M):

            sum_list = []

            for i in range(1, M):
                sum_list.append((math.sin(i * math.pi * currenttime) * eta[i]) / i * math.pi)

            currenttime += t
            first_ksi_arr.append(t * eta[0] + math.sqrt(2) * sum(sum_list))

        return first_ksi_arr

    def wiener2(eta, eta1, eta2):
        second_ksi_arr = []
        currenttime = 0

        for _ in range(M):

            sum_list = []

            for i in range(1, M):
                sum_list.append((eta1[i] * ((math.sin(2 * math.pi * i * currenttime)) / 
                (2 * math.pi * i))) + (eta2[i] * ((1 - math.cos(2 * math.pi * i * currenttime)) / (2 * math.pi * i))))
            currenttime += t
            second_ksi_arr.append(t * eta[0] + math.sqrt(2) * sum(sum_list))
        return second_ksi_arr

    #  графік реалізацій процесу

    plt.plot(np.arange(0, T, t), wiener1(eta))
    plt.show()
    plt.plot(np.arange(0, T, t), wiener2(eta, eta1, eta2))
    plt.show()

    
    # Перевірка точністі  для заданої кількості доданків у реалізації

    delta = 0.1
    epsilon = 0.05
    first_condition = pow(math.pi,(-2))*pow(delta,(-2))
    second_condition = math.exp(1 / 2) * delta * math.pi * math.sqrt(M) * math.exp((pow(delta,2) * pow(math.pi,2) * M)/(-2))

    if M > first_condition and second_condition <= epsilon:

        print("Model W1 approximates the process")

    else:

        print("Model W1 does not approximate the process")
    

    j1 = (mp.nsum(lambda k: k ** (-2) / (math.pi ** 2), [M + 1, math.inf]))
    j2 = math.sqrt(mp.nsum(lambda k: k ** (-4) / (math.pi ** 2), [M + 1, math.inf]))
    condition_teorema2 = pow((((delta - j1) / j2) + 1), 1/2) * math.exp((pow(delta,2)-j1) / ((-2) * j2))

    if M > first_condition and condition_teorema2 <= epsilon:

        print("Model W2 approximates the process")

    else:
        print("Model W2 does not approximate the process")


    # Оцінка середнього значення та дисперсії за 100 реалізацій венерівського процесу

    realization_mean = []
    realization_var = []
    realization = []

    # Побудова 100 реалізацій

    for _ in range(100):

        eta_ = np.random.normal(0, 0.1, M)
        
        realization.append(wiener1(eta_)) 

        # В кожній точці оцінити середнє значення та дисперсію
        realization_mean.append(mean(wiener1(eta_)))
        realization_var.append(variance(wiener1(eta_)))

    
    # Графік середнього значення і дисперсії

    plt.plot(np.linspace(0, 1, len(realization_mean)), realization_mean)
    plt.plot(np.linspace(0, 1, len(realization_mean)), realization_var)
    plt.legend(["Mean", "Variance"])
    plt.show()


    # Емпіричний закон розподілу ймовірностей часу першого виходу вінерівського процесу за заданий рівень
    a = 0.5
    time_lvl = []
    for i in realization:
          
        for j, value in enumerate(i):
            if value >= a:
                time_lvl.append(j * t)
                break

    sns.histplot(time_lvl, edgecolor = "black", color= "pink", bins=10)
    plt.show()






def start():
    print("Choose task you want to solve: \n1 - Poisson process\n2 - Wiener process")
    number = int(input(""))

    if number == 1:
        PuassonProcess()
        again()
    elif number == 2:
        WienerProcess()
        again()   
    else:
        print("Wrong number!")


def again():
    print("Would you like to try again?:\n1 - yes\n2 - no")
    choice = int(input(""))
    if choice == 1:
        start()
    
    elif choice == 2:
        return 0 


start()