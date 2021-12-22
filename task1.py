import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import derivative


def f(x):
    """
    Преобразовывает введенную функцию, которая в виде строки,
    к функции с привязанными математическими методами и парметром x.

    К примеру: f('sin(x)', 10) вернет функцию, которая
    высчитывает sin(10)

    :param x: str -- математическая функция в виде строки.

    :return: Выражение Python.

    """
    return eval(func,
                {'x': x,
                 'cos': math.cos, 'sin': math.sin, 'tg': math.tan,
                 'tan': math.tan, 'exp': math.exp, 'e': math.exp,
                 'log2': math.log2,
                 'ln': math.log, 'sqrt': math.sqrt})


def deriv(l_interval, r_interval, prec):
    """
    Дифференцирует функцию.

    :param l_interval: float -- левый конец интервала, на котором вы хотите
    продифференцировать функцию.
    :param r_interval: float -- правый конец интервала, на котором вы хотите
    продифференцировать функцию.
    :param prec: float -- шаг/точность, с которым вы хотите
    продифференцировать функцию.

    """
    n_1 = 0  # количество точек (int)

    deriv_x, deriv_y = [], []  # Значение точек х и
    # значение дифференцирования (list)

    for step_der in np.arange(l_interval, r_interval, prec):
        try:
            df = (f(step_der + prec) - f(step_der)) / prec  # дифференцирование
            # точек вручную (float)
        except ValueError:
            continue

        deriv_x.append(step_der)
        deriv_y.append(df)

        n_1 += 1

        try:
            df2 = derivative(f, step_der)  # дифференцирование встроенной
            # функцией scipy (float)
        except ValueError:
            print('{} | {} | {} | {} | {}'.format(n_1, round(step_der, 4),
                                                  df, None, None))
            continue

        print('{} | {} | {} | {} | {}'.format(n_1, round(step_der, 4), df,
                                              df2, df - df2))

    list_x, list_y = [], []  # список точек (list)

    for step_der in np.arange(l_interval, r_interval, 0.01):
        try:
            list_y.append(f(step_der))
        except ValueError:
            continue
        list_x.append(step_der)

    plt.title("Графики функции и производной f(x)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.scatter(list_x, list_y, color='b')
    plt.scatter(list_x, list_y, label="f(x)", color='b')
    plt.scatter(deriv_x, deriv_y, color='c')
    plt.scatter(deriv_x, deriv_y, label="f'(x)", color='c')
    plt.legend()

    plt.show()


def integ(l_interval, r_interval, count_n):
    """
    Интегрирует функцию.

    :param l_interval: float -- левый конец интервала, на котором вы хотите
    интегрировать функцию.
    :param r_interval: float -- правый конец интервала, на котором вы хотите
    интегрировать функцию.
    :param count_n: float -- шаг сетки.

    """
    s = 0  # площадь интегрируемой фигуры (float)
    value_h = (r_interval - l_interval) / count_n  # шаг сетки (float)
    for step_int in np.linspace(l_interval, r_interval, count_n + 1)[1:-1]:
        try:
            s += f(_)
            s += f(l_interval) / 2
            s += f(r_interval) / 2
        except ValueError:
            continue

    s *= value_h
    print('Значение интеграла функции на отрезке '
          '[{}; {}] равно: {}'.format(l_interval, r_interval, s))


trig_func_lst = ['cos', 'sin', 'tg', 'sec', 'cosec']

func = input('Введите функцию: ')
func = func.replace('^', '**')

operation = int(input('Выберите один из вариантов:\n'
                      '0 - дифференцировать функцию,\n'
                      '1 - интегрировать функцию. \n'))

for element in trig_func_lst:
    if func.lower().find(element) != -1:
        angle = int(input('Введите: \n'
                          '0 - если углы в радианах, \n'
                          '1 - если в градусах\n'))

l_int = float(input('Введите левую границу интервала: '))
r_int = float(input('Введите правую границу интервала: '))
if l_int > r_int:
    print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
    quit()

pr = 1
count = 1
if operation == 0:
    step_pr = int(input('Введите: \n'
                        '0, если вы хотите задать точность, \n'
                        '1, если шаг.\n'))
    if step_pr == 0:
        pr = float(input('Введите точность: '))
        deriv(l_int, r_int, pr)
    elif step_pr == 1:
        pr = float(input('Введите шаг: '))
        deriv(l_int, r_int, pr)
    else:
        print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
        quit()
elif operation == 1:
    count_step = int(input('Введите: \n'
                           '0 - если вы хотите задать количество отрезков, \n'
                           '1 - если шаг.\n'))
    if count_step == 0:
        count = int(input('Введите кол-во отрезков: '))
        integ(l_int, r_int, count)
    elif count_step == 1:
        count = round((r_int - l_int) /
                      float(input('Введите шаг: ')))
        integ(l_int, r_int, count)
    else:
        print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
        quit()
else:
    print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
    quit()
