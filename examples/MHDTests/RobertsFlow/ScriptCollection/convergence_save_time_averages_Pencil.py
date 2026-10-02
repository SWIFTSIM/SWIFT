
import numpy as np
import pandas as pd

results_directory_name = 'data_for_export_Pencil_one_mode'

def load_test_run_parameters():
    run_data = pd.read_csv("./"+results_directory_name + "/" + "run_parameters.csv", sep="\t")
    print(run_data)
    mask = run_data["status"] == "done"
    run_data = run_data[mask]
   # mask = run_data['run.in:kinematic_flow']=='Roberts-I'
   # run_data = run_data[mask]
   # mask = run_data['run.in:eta']<0.02
   # run_data = run_data[mask]
    return run_data

def find_growth_rate(time, B_field,time_interval):
    mask = ( ( time>time_interval[0] ) & ( time<time_interval[-1] ) )
    print(sum(mask))
    res, cov = np.polyfit(time[mask], np.log(B_field[mask]), 1,cov=True)
    return res, cov

def find_average(time, q, time_interval):
    mask = ( ( time>time_interval[0] ) & ( time<time_interval[-1] ) )
    mean_q = np.mean(q[mask])
    std_q = 2.0*np.std(q[mask])
    print('quantity:',q,'relative error:', std_q/mean_q)
    return mean_q, std_q

def load_dat_file(addr, time_interval):
    the_statistics = np.transpose(np.loadtxt(addr))
    
    # load data
    time = the_statistics[1] 
    urms = the_statistics[4] 
    brms = the_statistics[5] 
    jrms = the_statistics[6]
    abm = the_statistics[8]
    jbm = the_statistics[9]
    bmz = the_statistics[10]
    ujxbm = the_statistics[12]
    jbrms = the_statistics[13]
    jxbrms = the_statistics[14]
    j2b2m = the_statistics[15]

    # compute growth rate
    res,cov = find_growth_rate(time, brms,time_interval)
    lam_t = res[0]
    lam_t_error = 2.0*np.sqrt(cov[0][0])
    
    
    # compute mean quantities
    jrms_over_brms_t,jrms_over_brms_t_err = find_average(time, jrms/brms, time_interval)   
    jbm_over_brms2_t,jbm_over_brms2_t_err = find_average(time, jbm/brms**2, time_interval)   
    bmz_over_brms_t,bmz_over_brms_t_err = find_average(time, bmz/brms, time_interval)   
    ujxbm_over_brms2_t,ujxbm_over_brms2_t_err = find_average(time, ujxbm/(urms*brms**2), time_interval)   
    jbrms_over_brms2_t,jbrms_over_brms2_t_err = find_average(time, jbrms/brms**2, time_interval)   
    jxbrms_over_brms2_t, jxbrms_over_brms2_t_err = find_average(time, jxbrms/brms**2, time_interval)   
    j2b2m_over_brms2_t, j2b2m_over_brms2_t_err = find_average(time, j2b2m/brms**4, time_interval)   
    abm_over_brms2_t, abm_over_brms2_t_err = find_average(time, abm/brms**2, time_interval)

    keys_Pencil = ['lam','Bm/Brms','-WL/B2', 'Jrms/Brms', 'JB/B2', 'JBrms/B2']
    result = {'lam':lam_t,'lam_error':lam_t_error,'Bm/Brms': bmz_over_brms_t, '-WL/B2': -ujxbm_over_brms2_t, 'Jrms/Brms':jrms_over_brms_t, 'JB/B2':jbm_over_brms2_t,'JBrms/B2':jbrms_over_brms2_t, 'JB/JrmsBrms':jbm_over_brms2_t/jrms_over_brms_t,'<JB>rms/JrmsBrms':jbrms_over_brms2_t/jrms_over_brms_t,'<JxB>rms/JrmsBrms':jxbrms_over_brms2_t/jrms_over_brms_t, '<J2xB2>rms/(JrmsBrms)^2':j2b2m_over_brms2_t/(jrms_over_brms_t)**2, '-WL/(urms*brms*jrms)': -ujxbm_over_brms2_t/jrms_over_brms_t, '-<AB>/Brms2':-abm_over_brms2_t }

    return result

def process_info(run_data,time_interval):
    results = []
    for i in range(len(run_data)):
        run_data_slice = run_data.iloc[[i]]
        eta = run_data_slice["run.in:eta"].values[0]
        scheme = run_data_slice["Makefile.local:DERIV"].values[0]
        run_name = run_data_slice["run_name"].values[0]
        Nperp = int( run_name.split('_')[2])
        addr = (
            results_directory_name + '/'
            + str(run_data_slice["run_name"].values[0])
            + "/time_series.dat"
        )
        result = load_dat_file(addr, time_interval)
        result['eta'] = eta
        result['scheme'] = scheme
        result['run_name'] = run_name
        result['Nx'] = Nperp
        results.append(result)

    results_df = pd.DataFrame(results)
    return results_df

def save_results(data,addr):
    data.to_csv(addr, sep=" ", index=False,float_format="%.17g")
    return 


run_data = load_test_run_parameters()
results = process_info(run_data,time_interval=[97,100])

# save 2nd order Pencil
mask = results['scheme']=='deriv_2nd'
save_results(results[mask],'./data_for_export_Pencil_one_mode/5em2_2nd_rerun.txt')

# save 6th order Pencil
mask = results['scheme']=='deriv_6th'
save_results(results[mask],'./data_for_export_Pencil_one_mode/5em2_6th_rerun.txt')

# save 10th order Pencil
mask = results['scheme']=='deriv_10th'
save_results(results[mask],'./data_for_export_Pencil_one_mode/5em2_10th_rerun.txt')
