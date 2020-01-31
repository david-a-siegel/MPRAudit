def uncertainty(a,b,da,db):
    da = np.sqrt(da) #I'm inputting the variance, should be std dev.
    db = np.sqrt(db) #I'm inputting the variance, shoudl be std dev.
    return 1-a/b,((b*da)**2+(a*db)**2)/b**4

