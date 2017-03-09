__author__ = 'sebasanper'
# Eddy Viscosity wake model applied to horns rev.
import wake
from eddy_viscosity import ainslie
from math import sqrt, log, cos, sin
import time
from numpy import deg2rad

output = open('matrix_eddy.dat', 'w')
output.write('# This file has the wake deficit matrix per turbine per wind direction\n')
output2 = open('final_speed_eddy.dat', 'w')
output2.write('# This file has the deficit, wind speed and power at each turbine per wind direction.\n# Turbine number\tX-coordinate\tY-coordinate\tTotal speed deficit\tTotal wind speed\tWind direction angle\tPower produced\n')
layout = open('horns_rev.dat', 'r')
windrose = open('horns_rev_windrose2.dat', 'r')
draw = open('draw_horns_rev_eddy.dat', 'w')
draw.write('# This file has the turbines affected by the wake of one turbine at one direction.\n')
# draw2 = open('drawline.dat', 'w')
turb_data = open('turb17_data_eddy.dat', 'w')
direction = open('direction_efficiency_eddy.dat', 'w')
direction.write('# This file includes the efficiency of the whole farm by wind direction.\n# Wind direction angle\tFarm efficiency\n')

def determine_front(wind_angle, x_t1, y_t1, x_t2, y_t2):
    wind_angle = deg2rad(wind_angle)
    a = (x_t2 - x_t1) * cos(wind_angle) + (y_t2 - y_t1) * sin(wind_angle)
    if a > 0.0:
        return True, a
    else:
        return False, 0.0

def analysis():
    nt = 80  # Number of turbines
    D = 80.0  # Diameter
    layout_x = []
    layout_y = []
    for line in layout:
        columns = line.split()
        layout_x.append(float(columns[0]) / D)
        layout_y.append(float(columns[1]) / D)

    windrose_angle = []
    windrose_speed = []
    windrose_frequency = []
    for line in windrose:
        columns = line.split()
        windrose_angle.append(float(columns[0]))
        windrose_speed.append(float(columns[1]))
        windrose_frequency.append(float(columns[2]))

    layout.close()
    windrose.close()
    summation = 0.0

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 >= 4.0:
            return 19.7907842158 * U0 ** 2.0 - 74.9080669331 * U0 + 37.257967033  # From 4 to 11 m/s

    # def Cp(U0):
    #     return power(U0) / (0.5 * 1.225 * pi * r0**2.0 * U0**3.0)

    # for U0 in range(4, 20):
    #     summation = 0.0
    for wind in range(0, len(windrose_speed)):
    # for wind in range(0, 1):
        # U1 = windrose_speed[wind]  # Free stream wind speed
        # U0 = U1 * (70.0 / 10.0) ** 0.11  # Power or log law for wind shear profile
        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
        U0 = 8.5
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        for turbine in range(0, nt):
            parallel_distance = [0.0 for x in range(0, nt)]
            perpendicular_distance = [0.0 for x in range(0, nt)]
            flag = [False for x in range(nt)]
            for i in range(nt):
                flag[i], parallel_distance[i] = determine_front(angle3, layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i])
                perpendicular_distance[i] = wake.crosswind_distance(deg2rad(angle3), layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i])

            # Matrix with effect of each turbine <i = turbine> on every other turbine <j> of the farm
            for j in range(nt):
                if flag[j] is True:
                    if perpendicular_distance[j] <= 2.5:
                        wake_deficit_matrix[turbine][j] = ainslie(0.81, U0, parallel_distance[j], perpendicular_distance[j])
                    else:
                        wake_deficit_matrix[turbine][j] = 0.0
                elif flag[j] is False:
                    wake_deficit_matrix[turbine][j] = 0.0
                output.write('{0:f}\t'.format(wake_deficit_matrix[j][turbine]))
            output.write('\n')

        total_deficit = [0.0 for x in range(nt)]
        total_speed = [0.0 for x in range(nt)]
        for j in range(nt):
            for i in range(nt):
                total_deficit[j] += wake_deficit_matrix[i][j] ** 2.0
            total_deficit[j] = sqrt(total_deficit[j])
            total_speed[j] = U0 * (1.0 - total_deficit[j])
            output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(j, layout_x[j] * D, layout_y[j] * D, total_deficit[j], total_speed[j], int(angle), power(total_speed[j])))
        output2.write('\n')

        turb_data.write('{0:f}\t{1:f}\n'.format(angle, power(total_speed[14])))

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        efficiency = 0.0
        for l in range(nt):
            profit += power(total_speed[l])
            efficiency = profit * 100.0 / (nt * power(U0))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        # direction.write('{0:f}\t{1:f}\n'.format(angle, efficiency))
        summation += efficiency_proportion[wind]
    print 'total farm efficiency is {0:f} %'.format(summation)
        # print U0, summation

    turb_data.close()
    output.close()
    output2.close()
    draw.close()
    direction.close()

if __name__ == '__main__':
    start_time = time.time()
    analysis()
    print("--- %s seconds ---" % (time.time() - start_time))
