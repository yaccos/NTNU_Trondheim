from input_file import *
import classes
import datetime
now = datetime.datetime.now()


## Functions which prints data to file

def print_fozard(filename, delta_time):
    fi = open(filename, "w")


    print("MODEL USING FOZARD INPUT DATA",file=fi)
    print(now.strftime("%Y-%m-%d %H:%M"),file=fi)  
    print("Time used(min): %.1f" % delta_time, file=fi)

    print("\nINPUT",file=fi)
    print("%i timesteps of %.2f seconds to calculate %.1f hours" % (N, dt*60, final_time/60) ,file=fi)
    print("Bulk concentration: %.2f" % conc_bulk,file=fi)
    print("Box size: %i, %i, %i" % (Lx/vl, Ly/vl, Lz/vl) ,file=fi)
    print("Side length: %i" % vl,file=fi)
    print("Positive feedback Kq: %i" % Kq,file=fi)
    fi.close

def print_output(filename, biofilm):
    dat = open(filename, "w")
    # xyz, cs, cqsm, cqsi
    print("x y z mass/avg_mass_cell EPS subst qsm qsi",file=dat)
    for v in biofilm.vortex_arr:
       print(v.x, v.y, v.z, "%.0f" % (v.get_mass()/avg_mass_cell), v.eps_amount, v.conc_subst, v.conc_qsm, v.conc_qsi,file=dat)

    dat.close()


