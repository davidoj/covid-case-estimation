import pandas as pd
import numpy as np
from scipy.stats.distributions import gamma

ICL_STD = gamma(1/0.45**2,0,18.8*0.45**2)
ICL_ITS = gamma(1/0.86**2,0,5.1*0.86**2)

DEATHS_DAYS_S = np.array([ICL_STD.cdf(a+1)-ICL_STD.cdf(a) for a in range(100)])
S_DAYS = np.array([ICL_ITS.cdf(a+1)-ICL_ITS.cdf(a) for a in range(60)])

DEATHS_DAYS = np.convolve(DEATHS_DAYS_S,S_DAYS)


def normalize_jh_data(jh,name):
    jh['Country/Region'] = jh[['Country/Region','Province/State']].replace(np.nan,''""'').agg(' - '.join, axis=1).str.strip('- ')
    jh = jh.drop(columns=['Province/State','Lat','Long'])
    jh = jh.melt(id_vars=['Country/Region'],value_name=name)
    jh['variable'] = pd.to_datetime(jh['variable'])
    return jh

def get_jh_data():

	jhcc = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv')
	jhd = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv')


	jhcc, jhd = normalize_jh_data(jhcc,'Confirmed'), normalize_jh_data(jhd,'Deaths')

	out = jhd.set_index(['Country/Region','variable']).join(jhcc.set_index(['Country/Region','variable']))
	out = out.loc[~out.index.duplicated()]
	return out

# see heading "Align Confirmed Cases" below

def transformCases(s,shift=-7):
    shifted = s.diff().shift(shift)
    return shifted


# This isn't very accurate, and a better nowcast might be useful.
# However, wherever we use it we then take a convolutions with a distributions with 
# very small weight on the relevant days, and it's probably not orders of magnitude out
# which an exponential fit can be

def fillforward(orig,fill=7, pred=7):
    s = orig.copy()
    data_x = np.linspace(1,pred,pred)
    data_pred = np.linspace(pred+1,pred+1+fill,fill)
    try:
        s[-fill:] = np.poly1d(np.polyfit(data_x,
                                       s[-(pred+fill):-fill],1))(data_pred)
    except ValueError:
        print('Invalid data for linear fit', s[-(pred+fill):-fill])
        # In this case, we really don't know what cases are likely to do
        s[-fill:] = np.nan
        return s
    return s

# see heading "Comparison of Expected Deaths" below


def expectedDeaths(s,fatality=0.008,shift=7):
    s_pad = np.pad(s,(len(DEATHS_DAYS)-2,0),constant_values=0)
    cv = np.convolve(fillforward(s_pad,fill=shift),DEATHS_DAYS,'valid')
    pw = len(s)-len(cv)
    return fatality*np.pad(cv,(pw,0),"constant",constant_values=np.nan)

## see heading "Recovered, Infectious and Exposed Classes" below


def lik_r(i,mu=0.5):
    return np.exp(-mu*i)

norm = lik_r(np.arange(1,100)).sum()

def r(i):
    return lik_r(i)/norm

def R(ti):
    ti_pad = np.pad(ti,(40,0),'constant',constant_values=0)
    cv = np.convolve(ti_pad,r(np.arange(1,42)),'valid')
    pw = len(ti)-len(cv)
    return np.pad(cv,(pw,0),"constant",constant_values=0)

norm_I = lik_r(np.arange(1,100),0.2).sum()

def inf(i):
    return lik_r(i,0.2)/norm_I

def E2I(new_exposed):
    ne_pad = np.pad(new_exposed,(40,0),'constant',constant_values=0)
    cv = np.convolve(ne_pad,inf(np.arange(1,42)),'valid')
    pw = len(new_exposed)-len(cv)
    return np.pad(cv,(pw,0),"constant",constant_values=0)

# Calculate ascertainment, true infection rates, exposed and infectious classes and add as new columns
# Returns a much smaller dataframe, as 

def filtered_sum(m,indices):
    return m[indices].sum()

def ascertainment(csse_ds,fatality = 0.008,shift=7):
    csse_df = csse_ds.copy()
    csse_df['New confirmed shifted'] = csse_df['Confirmed'].groupby(level=0).transform(transformCases,shift=-1*shift)
    csse_df['New deaths'] = csse_df['Deaths'].groupby(level=0).transform(lambda x: x.diff())
    
    g = csse_df.groupby(level=0)
    
    csse_df['Expected deaths'] = g['New confirmed shifted'].transform(expectedDeaths,fatality=fatality,shift=shift)
    
    deaths = csse_df['New deaths']>=1
    last2 = csse_df.index.get_level_values(1)>= csse_df.index.get_level_values(1).max() - pd.Timedelta(days=14)
    
    indices = deaths #& last2
    
    csse_df['Ascertainment'] = np.nan
    
    csse_df.loc[indices,'Ascertainment'] = pd.to_numeric(csse_df.loc[indices,'Expected deaths']
                                            /csse_df.loc[indices,'New deaths'])
    
    asc_mean = (csse_df['Expected deaths'].groupby(level=0).transform(filtered_sum,indices)/
                csse_df['New deaths'].groupby(level=0).transform(filtered_sum,indices))
    csse_df['New cases true'] = (csse_df['New confirmed shifted']
                                             /asc_mean)

    g2 = csse_df.groupby(level=0)
    ever_exposed = g2['New cases true'].transform(np.cumsum)
    new_post_exposed= g2['New cases true'].transform(lambda x: E2I(x))
    post_exposed = new_post_exposed.groupby(level=0).transform(np.cumsum)
    csse_df['Exposed'] = ever_exposed - post_exposed
    csse_df['Recovered'] = new_post_exposed.groupby(level=0).transform(lambda x: np.cumsum(R(x)))
    csse_df['Infectious'] = post_exposed - csse_df['Recovered']
    return csse_df

def window_sum(a, n=3) :
    ret = a.cumsum()
    ret2 = a.shift(n).cumsum()
    return (ret-ret2)/n