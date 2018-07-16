from input_file_fozard import *
import classes
import datetime
now = datetime.datetime.now()

def print_fozard(biofilm, filename, delta_time):
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
    
    print("\nOUTPUT",file=fi)
    # xyz, cs, cqsm, cqsi
    for v in biofilm.vortex_arr:
       print(v.x, v.y, v.z, "%.0f" % v.get_mass(), v.eps_amount, v.conc_subst, v.conc_qsm, v.conc_qsi,file=fi)

    fi.close()


