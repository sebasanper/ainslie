__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Eddy Viscosity wake model applied to horns rev.
import wake
from eddy_viscosity_integrate import ainslie
from math import sqrt, log, cos, sin, tan
import time
from numpy import deg2rad

# output = open('matrix_eddy.dat', 'w')
# output.write('# This file has the wake deficit matrix per turbine per wind direction\n')
# output2 = open('final_speed_eddy.dat', 'w')
# output2.write('# This file has the deficit, wind speed and power at each turbine per wind direction.\n# Turbine number\tX-coordinate\tY-coordinate\tTotal speed deficit\tTotal wind speed\tWind direction angle\tPower produced\n')
layout = open('horns_rev.dat', 'r')
windrose = open('horns_rev_windrose2.dat', 'r')
# draw = open('draw_horns_rev_eddy.dat', 'w')
# draw.write('# This file has the turbines affected by the wake of one turbine at one direction.\n')
# draw2 = open('drawline.dat', 'w')
turb_data = open('turb14_ainslie.dat', 'w')
direction = open('direction_efficiency_ainslie.dat', 'w')
# direction.write('# This file includes the efficiency of the whole farm by wind direction.\n# Wind direction angle\tFarm efficiency\n')

def determine_front(wind_angle, x_t1, y_t1, x_t2, y_t2):
    wind_angle = deg2rad(wind_angle)
    a = (x_t2 - x_t1) * cos(wind_angle) + (y_t2 - y_t1) * sin(wind_angle)
    if a > 0.0:
        return a
    else:
        return 0.0

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

    def power(U):
        if U < 4.0:
            return 0.0
        elif U <= 25.0:
            return 0.0003234808 * U ** 7.0 - 0.0331940121 * U ** 6.0 + 1.3883148012 * U ** 5.0 - 30.3162345004 * U ** 4.0 + 367.6835557011 * U ** 3.0 - 2441.6860655008 * U ** 2.0 + 8345.6777042343 * U - 11352.9366182805
        else:
            return 0.0

    def Ct(U):
        return 0.0001923077 * U ** 4.0 + -0.0075407925 * U ** 3.0 + 0.096462704 * U ** 2.0 - 0.5012354312 * U + 1.7184749184


    def distance_to_front(x, y, theta, r):
        theta = deg2rad(theta)
        return abs(x + tan(theta) * y - r / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)

    # aver = [0.0 for x in range(nt)]
    for wind in range(0, len(windrose_angle)):
        # print wind
    # for wind in range(0, 1):
        U1 = windrose_speed[wind]  # Free stream wind speed
        U0 = U1 * (70.0 / 10.0) ** 0.11  # Power or log law for wind shear profile
        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
        # U0 = 8.5
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        distance = [[0.0 for x in range(2)] for x in range(nt)]
        total_deficit = [0.0 for x in range(nt)]
        total_speed = [U0 for x in range(nt)]

        for tur in range(nt):
            distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle, 100000000.0), tur]
        distance.sort()

        for turbine in range(nt):
            for num in range(turbine):
                total_deficit[distance[turbine][1]] += wake_deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0
            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
            total_speed[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
            parallel_distance = [0.0 for x in range(0, nt)]
            perpendicular_distance = [0.0 for x in range(0, nt)]
            for i in range(turbine + 1, nt):
                parallel_distance[distance[i][1]] = determine_front(angle3, layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]])
                perpendicular_distance[distance[i][1]] = wake.crosswind_distance(deg2rad(angle3), layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]])
                if perpendicular_distance[distance[i][1]] <= 1.7: ## 1.7 gives same results as a bigger distance, many times faster.
                    wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = ainslie(Ct(total_speed[distance[turbine][1]]), total_speed[distance[turbine][1]], parallel_distance[distance[i][1]], perpendicular_distance[distance[i][1]])
                else:
                    wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

        # for turbine in range(0, nt):
        #     parallel_distance = [0.0 for x in range(0, nt)]
        #     perpendicular_distance = [0.0 for x in range(0, nt)]
        #     flag = [False for x in range(nt)]
        #     for i in range(nt):
        #         flag[i], parallel_distance[i] = determine_front(angle3, layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i])
        #         perpendicular_distance[i] = wake.crosswind_distance(deg2rad(angle3), layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i])

            # Matrix with effect of each turbine <i = turbine> on every other turbine <j> of the farm

            #     output.write('{0:f}\t'.format(wake_deficit_matrix[j][turbine]))
            # output.write('\n')


        #     output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(j, layout_x[j] * D, layout_y[j] * D, total_deficit[j], total_speed[j], int(angle), power(total_speed[j])))
        # output2.write('\n')

        # for n in range(nt):
        #     aver[n] += power(total_speed[n]) / 360.0

        turb_data.write('{0:f}\n'.format(power(total_speed[14])))

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        for l in range(nt):
            profit += power(total_speed[l])
        efficiency = profit * 100.0 / (nt * power(total_speed[distance[0][1]]))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        direction.write('{0:f}\n'.format(profit))
        summation += efficiency_proportion[wind]
    print 'total farm efficiency is {0:f} %'.format(summation)
    # for n in range(nt):
        # turb_data.write('{0:f}\n'.format(aver[n]))
        # print U0, summation

    turb_data.close()
    # output.close()
    # output2.close()
    # draw.close()
    direction.close()

if __name__ == '__main__':
    start_time = time.time()
    analysis()
    print("--- %s seconds ---" % (time.time() - start_time))