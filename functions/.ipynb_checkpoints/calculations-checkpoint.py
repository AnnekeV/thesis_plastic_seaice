def round_down(x, decimal=1):
    '''Rounds down to specified number of decimals, default is one decimal'''
    return int(x*(10**decimal))/(10.**decimal)

def roll_average(data, windowsize):
    '''Smooths a data series by taking the average over specified window. note that it just works for xaray'''
    backupdata = data.copy()[:]
    data       = data.rolling(obs=windowsize, center=True, min_periods=1).mean()
    return data.where(~np.isnan(data), backupdata)
