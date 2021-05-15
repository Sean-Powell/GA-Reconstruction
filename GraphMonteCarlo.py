from matplotlib import pyplot as plt
from time import time
from os import mkdir
import numpy as np

COEFFICIENT_BIT_LENGTH = 5
EXPONENT_BIT_LENGTH = 3
FLOAT_PRECISION_LENGTH = 10
X_EXTENSION = 0.2
Y_UPPER_LIMIT = 300
Y_LOWER_LIMIT = 30
CL = 68 # Confidence limit in %

DATASET_PATH = "C:/Users/seanp/PycharmProjects/ThesisRewrite/Datasets/"
DATASET_NAME = "Hz_CC.txt"
MONTECARLO_PATH = "C:/Users/seanp/PycharmProjects/ThesisRewrite/monte_carlo_poly.txt"
x = []
y = []
s = []

time_str = str(time())
mkdir(time_str)

_float_bit_length = -1  # is automatically set by calculate_float_bit_length()
_term_bit_length = -1  #

def calculate_float_bit_length():
    global _float_bit_length
    number = '9' * FLOAT_PRECISION_LENGTH
    number = int(number)
    length = 0
    n = 2 ** length
    while n < number:
        length += 1
        n = 2 ** length
    _float_bit_length = length


def calculate_term_bit_length():
    global _term_bit_length
    length = 2  # for both signs
    length += (2 * _float_bit_length)
    length += COEFFICIENT_BIT_LENGTH
    length += EXPONENT_BIT_LENGTH
    _term_bit_length = length


def read_dataset():
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


def evaluate_chromosome(chromosome, value): # todo rewrite so that it only
    index = 0
    ans = 0
    while index < len(chromosome):
        coeff_sign = int(chromosome[index:index+1])
        index += 1
        coeff_i = chromosome[index: index + COEFFICIENT_BIT_LENGTH]
        index += COEFFICIENT_BIT_LENGTH
        coeff_d = chromosome[index: index + _float_bit_length]
        index += _float_bit_length

        try:
            exp_sign = int(chromosome[index:index + 1])
        except ValueError:
            exp_sign = 0
        index += 1
        exp_i = chromosome[index: index + EXPONENT_BIT_LENGTH]
        index += EXPONENT_BIT_LENGTH
        exp_d = chromosome[index: index + _float_bit_length]
        index += _float_bit_length

        if coeff_sign:
            coeff_sign = "-"
        else:
            coeff_sign = "+"

        if exp_sign:
            exp_sign = "-"  #
        else:
            exp_sign = "+"

        coeff = float(coeff_sign + binary_to_float(coeff_i, coeff_d))
        exp = float(exp_sign + binary_to_float(exp_i, exp_d))
        ans += (coeff * (value ** exp))

    return ans


def string_to_decimal(n):
    num = n
    dec_value = 0

    base1 = 1

    len1 = len(num)
    for i in range(len1 - 1, -1, -1):
        if num[i] == '1':
            dec_value += base1
        base1 = base1 * 2

    return dec_value


def binary_to_float(i, d):
    i = str(string_to_decimal(i))
    d = str(string_to_decimal(d))
    return i + "." + d


def plot_monte_carlo(chromosomes):
    # plt.plot(x, y, 'o')
    plt.errorbar(x, y, yerr=s, fmt='s', ecolor="tomato", elinewidth=0.5, capsize=1.5, ms=3, mew=0.5, markeredgecolor="tomato", color="white")
    func_x = np.linspace(0.001, max(x) + 0.2, 200)
    for c in chromosomes:
        func_y = []
        for i in func_x:
            func_y.append(evaluate_chromosome(c, i))

        plt.plot(func_x, func_y)

    plt.xlim(0, max(x) + X_EXTENSION)
    plt.ylim(Y_LOWER_LIMIT, Y_UPPER_LIMIT)
    plt.xlabel("Redshift")
    plt.ylabel("Distance modulus")
    plt.savefig(time_str + "/monte-carlo.png")
    plt.show()


def calculate_bands(chromosomes):
    alpha = CL / 100
    # plt.plot(x, y, 'o')
    plt.errorbar(x, y, yerr=s, fmt='s', ecolor="tomato", elinewidth=0.5, capsize=1.5, ms=3, mew=0.5, markeredgecolor="tomato", color="white")
    func_x = np.linspace(0.001, max(x) + 0.2, 200)
    bs_low, bs_up, bs_mu = [], [], []

    for i in func_x:
        y_bs_i = []
        for c in chromosomes:
            y_bs_i.append(evaluate_chromosome(c, i))

        p = (alpha + ((1.0 - alpha) / 2.0)) * 100
        upper = np.percentile(y_bs_i, p)

        p = ((1.0 - alpha) / 2.0) * 100
        lower = np.percentile(y_bs_i, p)

        bs_low.append(lower)
        bs_up.append(upper)
        bs_mu.append(np.percentile(y_bs_i, 50))

    plt.plot(func_x, bs_low, ls="--", c='blue')
    plt.plot(func_x, bs_up, ls='--', c='blue')
    plt.fill_between(func_x, bs_low, bs_up, alpha=0.6, edgecolor='royalblue',
                     facecolor='royalblue', linewidth=0.5, zorder=1, label='1 $\sigma$')
    plt.plot(func_x, bs_mu, 'k:', label='Mean Fit', c='black')
    plt.legend(loc='best')
    plt.xlim(0, max(x) + X_EXTENSION)
    plt.ylim(Y_LOWER_LIMIT, Y_UPPER_LIMIT)
    plt.xlabel("Redshift")
    plt.ylabel("Distance modulus")
    plt.savefig(time_str + "/monte-carlo-bands.png")
    plt.show()


def chromsome_to_string(chromosome, value):
    index = 0
    ans = ""
    while index < len(chromosome):
        coeff_sign = int(chromosome[index:index + 1])
        index += 1
        coeff_i = chromosome[index: index + COEFFICIENT_BIT_LENGTH]
        index += COEFFICIENT_BIT_LENGTH
        coeff_d = chromosome[index: index + _float_bit_length]
        index += _float_bit_length

        try:
            exp_sign = int(chromosome[index:index + 1])
        except ValueError:
            exp_sign = 0
        index += 1
        exp_i = chromosome[index: index + EXPONENT_BIT_LENGTH]
        index += EXPONENT_BIT_LENGTH
        exp_d = chromosome[index: index + _float_bit_length]
        index += _float_bit_length

        if coeff_sign:
            coeff_sign = "-"
        else:
            coeff_sign = "+"

        if exp_sign:
            exp_sign = "-"
        else:
            exp_sign = "+"

        coeff = float(coeff_sign + binary_to_float(coeff_i, coeff_d))
        exp = float(exp_sign + binary_to_float(exp_i, exp_d))

        if coeff_sign == "+":
            ans = ans + "+" + str(coeff) + "x^" + str(exp)
        else:
            ans = ans + "-" + str(coeff) + "x^" + str(exp)
        # ans += (coeff * (value ** exp))

    return ans


def plot_chrom(chrom):
    plt.plot(x, y, 'o')
    plt.errorbar(x, y, yerr=s, fmt='s', ecolor="tomato", elinewidth=0.5, capsize=1.5, ms=3, mew=0.5, markeredgecolor="tomato", color="white")
    func_x = np.linspace(min(x), max(x) + 0.2, 200)

    func_y = []
    for i in func_x:
        func_y.append(evaluate_chromosome(chrom, i))

    plt.plot(func_x, func_y)

    plt.xlim(0, max(x) + X_EXTENSION)
    plt.ylim(Y_LOWER_LIMIT, Y_UPPER_LIMIT)
    plt.xlabel("Redshift")
    plt.ylabel("Distance modulus")
    plt.title("Monte-Carlo Simulation " + DATASET_NAME.split(".txt")[0])
    plt.savefig(time_str + "/monte-carlo-chromosome.png")
    plt.show()


read_dataset()
calculate_term_bit_length()
calculate_float_bit_length()

file_f = open(MONTECARLO_PATH)
chroms = []

lowest_chi = 9999999999999999999999999999999999
lowest_chrom = -1

for line in file_f:
    parts = line.split(",")
    chroms.append(parts[0])
    if float(parts[1]) < lowest_chi:
        lowest_chi = float(parts[1])
        lowest_chrom = parts[0]


plot_monte_carlo(chroms)
calculate_bands(chroms)
f = open(time_str + "/best_chrom.txt", "w")
f.write(chromsome_to_string(lowest_chrom, x))
print(chromsome_to_string(lowest_chrom, x))
f.close()
