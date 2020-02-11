import numpy as np

def naive_data_cleaning(f):
    f = np.array(f)
    m_all = len(f)
    mu = np.mean(f)
    sigma = np.std(f)
    f_min = np.min(f)
    f_max = np.max(f)
    z = 3.
    xl = mu - z*sigma
    xh = mu + z*sigma

    f_cleaned = f[f>xl]
    m_low = m_all - len(f_cleaned)
    f_cleaned = f_cleaned[f_cleaned<xh]
    m_high = m_all - m_low - len(f_cleaned)

    mu = np.mean(f_cleaned)
    sigma = np.std(f_cleaned)
    r_low = float(m_low) / float(m_all)
    r_high = float(m_high) / float(m_all)


    return f_cleaned, mu, sigma, r_low, r_high, [f_min, f_max]

def scott_rule(f, sigma=None):
    if sigma == None:
        sigma = np.std(f)
    N = len(f)
    width = 3.49*sigma*N**(-1./3.) # Scott's rule
    bins_num = np.ceil((np.max(f) - np.min(f))/width).astype(int)
    return bins_num

def sturges_rule(f):
    return np.ceil(1+np.log2(len(f))).astype(int)

def bins_rule(f, name):
    if name=='scott':
        bins_num = scott_rule(f)
    elif name == 'sturges':
        bins_num = sturges_rule(f)
    else:
        bins_num = int(name)
    return bins_num

def get_sample_points(f, N):
    f_min = np.min(f)
    f_max = np.max(f)
    x = np.linspace(f_min, f_max, N)
    dx = x[1] - x[0]
    return x, dx

def bootstrap_mean(x, reps=599, cl=0.95):
    x = x.ravel()
    n = len(x)
    alpha = (1.-cl)/2.
    mean_list = list()
    for i in range(reps):
        xb = np.random.choice(x, n)
        mean_list.append(np.mean(xb))
    mean = np.mean(mean_list)
    lb = np.percentile(mean_list, alpha*100)
    ub = np.percentile(mean_list, (1-alpha)*100)
    result = {'mean':mean, 'CI':[lb,ub], 'data':mean_list}


    return result