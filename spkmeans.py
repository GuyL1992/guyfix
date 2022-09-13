import pandas as pd
import numpy as np
import sys
import spkmeans
from enum import Enum

class Goals(Enum):
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"
    SPK = "spk"

def initial_points(file_name):

    """ call this method to read the input file and translate it to nested list represents the observations 
    Args:
        file_name (string): input file name (path)

    Returns:
        list[list[float]], int, int: 3 returned values, one represents the observed vectors (nested list) and two int values
        to represent the dimentions (n & d)
    """

    with open(file_name,'r') as file:

        def parse_row(row):
            return list(map(lambda point: float(point),row.split(',')))
        observations = list(map(lambda row: parse_row(row), file))

    return observations, len(observations), len(observations[0])

### K_means_pp process and the required help functions ###

def calc_distance(vec1, vec2): # calculate distance between two vectors 
    return sum((vec1 - vec2)**2)
          
def kmeans_pp(k,data_arr,n,d):
    """_summary_

    Args:
        k (int): number of clusters required 
        data_arr (list[list[double]]): input observations 
        n (int): first dimention
        d (int): second dimention

    Returns:
        list[int], dataFrame: list of initiate k index, and the actual init k vectors as df
    """
    np.random.seed(0)
    first_centroid = np.random.choice(n)
    centroids_index = []
    p_list =  np.zeros((n), dtype=np.float64) 
    d_list =  np.full(n, np.inf) 
    centroids_index.append(first_centroid)
    d_sum = 0
    j = 1
    while j < k:
        new_centroid = data_arr.loc[centroids_index[j-1]]
        for i in range(n):
            curr_distance = calc_distance(data_arr.loc[i], new_centroid)
            if curr_distance < d_list[i]:
                d_list[i] = curr_distance
        d_sum = d_list.sum()

        for i in range(n):
            p_list[i] = d_list[i] / d_sum 
        next_centroid = np.random.choice(n, size=None, p=p_list)
        centroids_index.append(next_centroid)
        j += 1
    
    return centroids_index, data_arr.loc[centroids_index]


def recieve_input(): 
    """ recieve and validate the user input. assigned it to the appropriate variables

    Returns:
        int, str, str: k desired by the user - might not be provided, the desired flow, the input file path
    """
    arguments_size = len(sys.argv)
    
    if arguments_size == 4:
        k_float = float(sys.argv[1])
        flow = sys.argv[2]
        input = sys.argv[3]      

    else :
        print("Invalid Input!")
        exit()
    
    return k_float, flow, input

def invalid_input(): 
    print("Invalid Input!")
    return True   

if __name__ == "__main__":

    k_float,goal,input_path = recieve_input()

    try:
        observations, n, d = initial_points(input_path)
    except:
        invalid_input()

    message = False
    try:
        if goal ==  (Goals.SPK.value):
            k = int(k_float)   
            if k_float != k or k < 0:
                message = invalid_input()
                exit()

            new_data = spkmeans.get_new_data(n,d,k,observations)
            prev_k = k
            k = len(new_data[0])
            if prev_k >= len(new_data):
                message = invalid_input()
                exit()

            df = pd.DataFrame(new_data)
            init_index, init_centroids = kmeans_pp(k,df,n,k)
            any(map(lambda i: print(init_index[i], end=",") if i < len(init_index) - 1 else print(init_index[i]),range(len(init_index))))
            init_centroids = init_centroids.values.tolist()
            spkmeans.kmeans(n,k,k,new_data,init_centroids)

        elif goal == (Goals.WAM.value):
            spkmeans.wam(n,d,observations)

        elif goal == (Goals.DDG.value):
            spkmeans.ddg(n,d,observations)

        elif goal == (Goals.LNORM.value):
            spkmeans.lnorm(n,d,observations)

        elif goal == (Goals.JACOBI.value):
            spkmeans.jacobi(n,d,observations)
        else:
            message = invalid_input()
            exit()
    except:
        print("An Error Has Occured") if not message else None
        exit()

