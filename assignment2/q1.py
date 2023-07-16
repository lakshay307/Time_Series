function[X, Y, Z, C, PDOP] = GNSS_A1_5(
    PR, sat_coord_1, sat_coord_2, sat_coord_3, sat_coord_4, sat_coord_5, approx_station)

# % Check if optional inputs are specified

X = approx_station(1)
Y = approx_station(2)
Z = approx_station(3)

P = eye(5, 5)

for i = 1:15

   # % inital computed peusorange
    ini_pr_sat1 = sqrt((sat_coord_1(1)-X) ^ 2 + (sat_coord_1(2)-Y)
                       ^ 2 + (sat_coord_1(3)-Z) ^ 2) + sat_coord_1(4)
    ini_pr_sat2 = sqrt((sat_coord_2(1)-X) ^ 2 + (sat_coord_2(2)-Y)
                       ^ 2 + (sat_coord_2(3)-Z) ^ 2) + sat_coord_2(4)
    ini_pr_sat3 = sqrt((sat_coord_3(1)-X) ^ 2 + (sat_coord_3(2)-Y)
                       ^ 2 + (sat_coord_3(3)-Z) ^ 2) + sat_coord_3(4)
    ini_pr_sat4 = sqrt((sat_coord_4(1)-X) ^ 2 + (sat_coord_4(2)-Y)
                       ^ 2 + (sat_coord_4(3)-Z) ^ 2) + sat_coord_4(4)
    ini_pr_sat5 = sqrt((sat_coord_5(1)-X) ^ 2 + (sat_coord_5(2)-Y)
                       ^ 2 + (sat_coord_5(3)-Z) ^ 2) + sat_coord_5(4)

#    % misclosure matrix
    k = PR - [ini_pr_sat1
              ini_pr_sat2
              ini_pr_sat3
              ini_pr_sat4
              ini_pr_sat5]
    # % misclosure

    first_row = [-(sat_coord_1(1) - X)/ini_pr_sat1, -(sat_coord_1(2) -
                                                      Y)/ini_pr_sat1, -(sat_coord_1(3) - Z)/ini_pr_sat1, 1]

    second_row = [-(sat_coord_2(1) - X)/ini_pr_sat2, -(sat_coord_2(2) -
                                                       Y)/ini_pr_sat2, -(sat_coord_2(3) - Z)/ini_pr_sat2, 1]

    third_row = [-(sat_coord_3(1) - X)/ini_pr_sat3, -(sat_coord_3(2) -
                                                      Y)/ini_pr_sat3, -(sat_coord_3(3) - Z)/ini_pr_sat3, 1]

    fourth_row = [-(sat_coord_4(1) - X)/ini_pr_sat4, -(sat_coord_4(2) -
                                                       Y)/ini_pr_sat4, -(sat_coord_4(3) - Z)/ini_pr_sat4, 1]

    fifth_row = [-(sat_coord_5(1) - X)/ini_pr_sat5, -(sat_coord_5(2) -
                                                      Y)/ini_pr_sat5, -(sat_coord_5(3) - Z)/ini_pr_sat5, 1]

#   %  jacobian matrix
    J = [first_row
         second_row
         third_row
         fourth_row
         fifth_row]
    % design matrix

    N = J'* P*J
    % normal equation

    U = J'* P*k

    X_mat = inv(N)*U

    X = X + X_mat(1)
    Y = Y + X_mat(2)
    Z = Z + X_mat(3)
    C = X_mat(4)
# end

a_priori = 1
cof_par = inv(N)
% cofactor matrix of parameters
covar_par = a_priori * cof_par
% variance covariance of parameters

# to compute the PDOP

PDOP = sqrt(trace(covar_par(1: 3, 1: 3)))

# end
