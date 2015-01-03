import pickle
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot, animation


if __name__ == "__main__":
    data = pickle.load(open("multi_gmc_timings.dat", "rb"))
    number_of_clouds = [i[0] for i in data]
    time_deltas = [i[1] for i in data]

    fig = pyplot.figure(figsize=[10,10])

    pyplot.scatter(number_of_clouds, time_deltas, c='black')
    pyplot.xlabel("number of gas clouds")
    pyplot.ylabel("seconds")
    pyplot.title("computation time of 4Myr simulations")

    pyplot.savefig('multi_gmc_timings.png')
    pyplot.close()
                

