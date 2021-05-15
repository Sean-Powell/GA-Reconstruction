import tensorflow as tf
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np

DATASET_PATH = "C:/Users/seanp/PycharmProjects/ThesisRewrite/Datasets/"
DATASET_NAME = "Hz_CC.txt"
MODEL_PATH = "C:/Users/seanp/PycharmProjects/ThesisRewrite/two_exp_single_relu/"

x = []
y = []
s = []

data_set_names = ["Hz_CC.txt", "Hz_CC_BAO.txt", "Pantheon_calibration_HW.txt", "salt_light_curve_constitution.txt",  "SCPUnion_mu_vs_z.txt", "SCPUnion2.1_mu_vs_z.txt"]
table_names = ["Hz-CC", "Hz-CC-BAO", "Pantheon-Calibration", "SALT-Constitution", "SCP-Union", "SCP-Union-2.1"]
x_values = []
s_values = []

def read_dataset():
    global x, y, s
    x = []
    y = []
    s = []

    f = open(DATASET_PATH + DATASET_NAME)
    i = 0

    for line in f:
        line = line.strip()
        if DATASET_NAME == "Hz_CC.txt" or DATASET_NAME == "Hz_CC_BAO.txt" or DATASET_NAME == "Pantheon_calibration_HW.txt":
            line_split = line.split("\t")
            x.append(float(line_split[0]))
            y.append(float(line_split[1]))
            s.append(float(line_split[2]))
        elif DATASET_NAME == "salt_light_curve_constitution.txt":
            if i < 1:
                i += 1
            else:
                line_split = line.split(" ")
                x_string = line_split[1]
                index = x_string.find("(")
                x_value = x_string[:index]
                x.append(float(x_value))

                y_string = line_split[5]
                index = x_string.find("(")
                y_value = y_string[:index + 1]
                s_value = y_string[index + 2: len(y_string) - 1]
                y.append(float(y_value))
                s.append(float(s_value))
        elif DATASET_NAME == "SCPUnion_mu_vs_z.txt":
            if i < 4:
                i += 1
            else:
                line_split = line.split("\t")
                x.append(float(line_split[1]))
                y.append(float(line_split[2]))
                s.append(float(line_split[3]))
        elif DATASET_NAME == "SCPUnion2.1_mu_vs_z.txt":
            if i < 5:
                i += 1
            else:
                line_split = line.split("\t")
                x.append(float(line_split[1]))
                y.append(float(line_split[2]))
                s.append(float(line_split[3]))


print(DATASET_NAME)

hist_list = []

for data_set in data_set_names:
    DATASET_NAME = data_set
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(units=1, activation='linear', input_shape=[1]))
    model.add(tf.keras.layers.Dense(units=1024, activation='exponential'))
    model.add(tf.keras.layers.Dense(units=1024, activation='exponential'))
    model.add(tf.keras.layers.Dense(units=1024, activation='relu'))
    model.add(tf.keras.layers.Dense(units=2, activation='linear'))
    model.compile(loss='mse', optimizer='adam')

    model.summary()
    read_dataset()
    x = np.array(x)
    y = np.array(y)
    s = np.array(s)
    data = np.array(list((zip(y, s))))
    hist = model.fit(x, data, epochs=1000, batch_size=16, verbose=1)
    hist_list.append(hist)
    fig, ax = plt.subplots()
    # plt.scatter(x[::1], y[::1])
    plt.errorbar(x[::1], y[::1], yerr=s[::1], fmt='s', ecolor="tomato", elinewidth=0.5, capsize=1.5, ms=3, mew=0.5,
                 markeredgecolor="tomato", color="white", zorder=2)
    func_x = np.linspace(0.001, max(x) + 0.05, 10000)
    predicted_values = model.predict(func_x)

    predicted_values_zero = model.predict(np.array([0]))
    print(predicted_values_zero)
    x_values.append(predicted_values_zero[0][0])
    s_values.append(predicted_values_zero[0][1])

    y_prediction = predicted_values[:, 0]
    s_prediction = predicted_values[:, 1]

    plt.fill_between(func_x, y_prediction - 2 * s_prediction, y_prediction + 2 * s_prediction, alpha=0.2,
                     edgecolor='royalblue', facecolor='royalblue', linewidth=0.5, zorder=1)
    plt.fill_between(func_x, y_prediction - s_prediction, y_prediction + s_prediction, alpha=0.6, edgecolor='royalblue',
                     facecolor='royalblue', linewidth=0.5, zorder=1)

    plt.plot(func_x, y_prediction, 'blue', zorder=3)
    plt.xlabel("Redshift")
    plt.ylabel("Distance Modulus")
    plt.xlim(left=0, right=(max(x) + 0.05))
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    plt.savefig(MODEL_PATH + DATASET_NAME.split(".txt")[0] + "_ANN.png")
    plt.show()

    model.save(MODEL_PATH + DATASET_NAME.split(".txt")[0])
    model.save_weights(MODEL_PATH + DATASET_NAME.split(".txt")[0] + "_Weights")

plt.figure()

for i in range(len(hist_list)):
    hist = hist_list[i]
    loss = hist.history['loss']
    epochs = range(1, len(loss) + 1)
    plt.plot(epochs, loss, label=data_set_names[i].split(".txt")[0])

plt.ylabel('loss')
plt.xlabel('epoch')
plt.yscale('log')
plt.legend(loc=7)
plt.xlim(left=0, right=1000)
plt.savefig(MODEL_PATH + "loss_graph.png")
plt.show()
plt.close()

f = open(MODEL_PATH + "x_values.txt", "w")
f.write("\hline\n")
f.write("Data set names & Y value at Z=0 & \(\sigma\)\\\\\n")
f.write("\hline\n")
for i in range(len(x_values)):
    x_value = x_values[i]
    s_value = s_values[i]
    name= table_names[i]
    f.write(str(name) + " & " + str(x_value) + " & " + str(s_value) + "\\\\\n")
f.write("\hline\n")
f.close()

