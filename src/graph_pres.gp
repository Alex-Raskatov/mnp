count = 1000
step = 10
max_z = 0.05

set xrange [0:100]
set yrange [0:100]
set cbrange [0:max_z]

set size ratio -1
set view map

set term gif animate delay 10 size 2000, 1000

set output "graph/part3_kin.gif"

do for [n = 0 : count - 1 : step] {
    filename = sprintf("data_kin/out_%03d.dat", n)
    splot filename with image title ""
}