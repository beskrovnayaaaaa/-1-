import scipy
import pywt
import os
import random
import math
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pywt import wavedec
from scipy import integrate
from scipy.fft import fft, fftfreq
from numpy import arange
from pylab import *
from numpy import *
from scipy import *
from sympy import *


def file(file_name):
    """
    Преобразует введенный файл, состоящий из двух столбцов x и y,
    в список списков из xi,yi.
    :param file_name: csv - csv файл
    :return: список списков из xi,yi
    """
    points_file = pd.read_csv(file_name, names=['x', 'y'])
    points_1 = points_file.values.tolist()
    return points_1


def c_n(func, st_1, lv):
    """
    Нахождение коээфициента несинусоидальности
    :param func: list - список значений аппроксимации или детализации
    :param st_1: str - метод нахождения, например 'db25' или 'haar'
    :param lv: int - уровень вейвлета
    :return: коэффициент несинусоидальности
    """
    a = []
    k = 0
    k1 = 0
    cf_n = []
    for i in range(1, lv):
        cf_n = []
        lv = i
        cf_n = wavedec(func, st_1, lv)
        a.append(np.std(cf_n[1]))
        k += a[i-1]**2
    k1 = math.sqrt(k) / np.std(cf_n[0])
    return k1


def cm_to_inch(value):
    """
    Преобразует дюймы в сантиметры
    :param value: int - дюймы
    :return: float - см
    """
    return value/2.54


def start(x_1, y_1):
    """
    Построение начального графика по файлу csv.
    :param x_1: list - список x
    :param y_1: list - список y
    :return: строит график
    """
    plt.figure(figsize=(cm_to_inch(40), cm_to_inch(30)))
    subplot(3, 1, 1)
    plt.plot(x_1, y_1, 'pink', linewidth=2, label='Изначальная функция')
    grid()
    legend(loc='best')


def daub(x, y, level):
    """
    Функция находит коэффициенты аппроксимации или детализации,
    используя метод Daubechies
    :param x: list - список x
    :param y: list - список y
    :param level: int - уровень вейвлета
    :return: строит графики для коэффициентов аппроксимации и деталиции заданного уровня
    """
    st = 'db25'
    coeffs = wavedec(y, st, level=level)
    print('level=', level+1)
    plt.figure(figsize=(cm_to_inch(40), cm_to_inch(30)))
    subplot(3, 1, 1)
    plt.plot(coeffs[0], 'b', linewidth=2,
             label='Коэффициенты аппроксимации нижних частот')
    grid()
    legend(loc='best')
    if (level == 0):
        print("График детализации нельзя вывести")
    else:
        coeffs1 = wavedec(y, st, level=level)
        print('level=', level+1)
        subplot(3, 1, 2)
        plt.plot(coeffs1[1], 'r', linewidth=2,
                 label='Коэффициенты детилизации высоких частот')
        grid()
        legend(loc='best')
        show()
    return ''


def haar(x, y, level1):
    """
    Функция находит коэффициенты аппроксимации или детализации,
    используя метод Haar
    :param x: list - список x
    :param y: list - список y
    :param level1: int - уровень вейвлета
    :return: строит графики для коэффициентов аппроксимации и деталиции заданного уровня
    """
    st1 = 'haar'
    coeffs1 = wavedec(y, st1, level=level1)
    print('level=', level1+1)
    plt.figure(figsize=(cm_to_inch(40), cm_to_inch(30)))
    subplot(3, 1, 1)
    plt.plot(coeffs1[0], 'b', linewidth=2,
             label='Коэффициенты аппроксимации нижних частот')
    grid()
    legend(loc='best')
    if (level1 == 0):
        print("График детализации нельзя вывести")
    else:
        coeffs2 = wavedec(y, st, level=level1)
        print('level=', level1+1)
        subplot(3, 1, 2)
        plt.plot(coeffs2[1], 'r', linewidth=2,
                 label='Коэффициенты детилизации высоких частот')
        grid()
        legend(loc='best')
        show()
    return ''


def s(c, y):
    """
    :param c: list - значения либо аппроксимации, либо детализации
    :param y: list - список изначальных y
    :return: значение квадратичного отклонения
    """
    ss = 0
    for i in range(0, len(c)):
        ss += (c[i] - y[i])**2
    return ss**0.5/len(c)


def wav(y, x, name):
    """
    :param y: list - список изначальных x
    :param x: list - список изначальных y
    :param name: метод нахождения вейвлетов ('mexh'- на основе Мексиканской шляпы,
    gaus1 - на основе Гаусса)
    :return: графики 4 уровней вейвлетов и коэффициент несинусоидальности
    """
    print('Вейвлет', name)
    wav = pywt.ContinuousWavelet(name)
    coef, frec = pywt.cwt(y, arange(1, 30), wav)
    d1 = s(coef[0], y)
    plt.plot(coef[0])
    plt.title('level=1')
    plt.show()
    coef, frec = pywt.cwt(coef[0], arange(1, 30), wav)
    d2 = s(coef[0], y)
    plt.plot(coef[0])
    plt.title('level=2')
    plt.show()
    coef, frec = pywt.cwt(coef[0], arange(1, 30), wav)
    d3 = s(coef[0], y)
    plt.plot(coef[0])
    plt.title('level=3')
    plt.show()
    coef, frec = pywt.cwt(coef[0], arange(1, 30), wav)
    d4 = s(coef[0], y)
    a4 = s(coef[1], y)
    plt.plot(coef[0])
    plt.title('level=4')
    plt.show()
    print((d1**2+d2**2+d3**2+d4**2)**0.5/a4)


def FFT(x, y):
    """
    Функция строит графики спектограммы, изначального графика и восстановленной
    функции, считает дисперсию и коэффициент несинусоидальности сигнала
    методом быстрого Фурье преобразования
    :param x: list - список изначальных x
    :param y: list - список изначальных y
    :return: ''
    """
    yf = rfft(y)
    plt.title('Спектрограмма')
    xf = rfftfreq(len(y), x[-1]/len(y))
    y_new = irfft(yf)
    plt.plot(xf, abs(yf)/len(y))
    plt.xlim([0, 10])
    plt.show()
    plt.title('Изначальный график')
    plt.plot(x, y)
    plt.show()
    plt.title('Восстановленная функция')
    plt.plot(x, y_new)
    plt.show()
    d = 0
    for i in range(len(yf)):
        d += y[i]-y_new[i]
    d = d/len(y)
    print('Дисперсия: ', d)
    S = scipy.integrate.trapz(abs(yf)/len(y), xf)
    x = xf
    y = list((abs(yf)/len(y)))
    max_y = max(y)
    index_max_y = y.index(max(y))
    if index_max_y == 0:
        a = max_y
        a_i = 0
    else:
        for i in range(index_max_y-1, -1, -1):
            if (y[i] < y[i-1]) & (y[i] < y[i+1]):
                a = y[i]
                a_i = i
                break
    for i in range(index_max_y+1, len(y)-1):
        if (y[i] < y[i-1]) & (y[i] < y[i+1]):
            b = y[i]
            b_i = i
            break
    S0 = scipy.integrate.trapz([y[i] for i in range(a_i, b_i+1)],
                         [x[i] for i in range(a_i, b_i+1)])
    Si = S - S0
    K = Si / S0
    print ("Коэффициент несинусоидальности сигнала:", K)
    return ''


name = input("Укажите путь к файлу в формате: \
    '/Users/Vlada/Downloads/имя_файла.csv' - ")
points = file(name)
x = [points[i][0] for i in range(len(points))]
print(len(x))
y = [points[i][1] for i in range(len(points))]


start(x, y)
print("На основе Daubechies")
level = 4
for i in range(0, level):
    daub(x, y, i)
st = 'db25'
print("Коэффициент несинусоидальности для вейвлет-анализа")
print(10*c_n(y, st, 4))
print("На основе Haar")
level1 = 4
for i in range(0, level1):
    haar(x, y, i)
st1 = 'haar'
print("Коэффициент несинусоидальности для вейвлет-анализа")
print(c_n(y, st1, 4))
print("На основе Мексиканская шляпа")
wav(y, x, 'mexh')
print("На основе Гаусса")
wav(y, x, 'gaus1')
print("Быстрое преобразование Фурье")
FFT(x, y)
