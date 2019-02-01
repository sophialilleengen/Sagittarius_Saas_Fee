def cosmo_cfunc(M200,h):
    #From Dutton & Maccio 2014:
    c = 10.**(0.905 - 0.101 * (np.log10(M200*h)-12.))
    return c

def corenfw_mass(r,M200,c,oden,tSF,Rhalf,rhocrit,eta,kappa):
    #Assumes input arrays in Msun, kpc:
    G = 6.67e-11
    kpc = 3.086e19
    Msun = 1.989e30
    Gyr = 365.*24.*60.*60.*1e9

    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c
    rhos=rhocrit*deltachar
    rhoanal = rhos/((r/rs)*(1.+(r/rs))**2.)
    manal = M200 * gcon * (np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))

    print('R_200:', rv)
    print('coreNFW scale length:', rs)
    print('coreNFW density parameter:', rhos)

    Mrs = M200 * gcon * (np.log(2.0)-0.5)
    tdyn_rs = \
        2.*np.pi*np.sqrt((rs*kpc)**3./G/(Mrs*Msun))/Gyr

    rc = eta * Rhalf

    x = r/rc
    f = (np.exp(x) - np.exp(-x))/(np.exp(x)+np.exp(-x))
    xx = tSF/tdyn_rs * kappa
    if (tSF > 0.0):
        n = (np.exp(xx) - np.exp(-xx))/(np.exp(xx)+np.exp(-xx))
    else:
        n = -tSF
    my_manal = manal*f**n

    my_rhoanal = rhoanal*f**n + \
        1.0/(4.*np.pi*r**2.*rc)*manal*(1.0-f**2.)*n*f**(n-1.0)

    return my_manal

#Simple program to calculate the properties of the NFW profile
#given an input M200. c200 is calculated assuming the mean
#c200-M200 relation in LCDM taken from Dutton & Maccio 2014.
if __name__ == "__main__":
    #Import plots library:
    import numpy as np
    import pylab as plt
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    #Parameters:
    h = 0.7          #Hubble parameter
    oden = 200.0     #Overdensity parameter (200=>M200 etc)
    tSF = 0.01       #Makes everything NFW (no core)
    Rhalf = 1.0      #Sets the core size, but we don't have one => ignore
    rhocrit = 135.05 #Critical density of Universe in Msun/kpc^3
    eta = 1.75       #Furthre core parameters => don't need
    kappa = 0.04     #

    #M200 in Msun
    M200 = 6.0e10

    #Calculate c200 using Dutton & Maccio:
    c200 = cosmo_cfunc(M200,h)

    #Calculate analytic NFW profile with the above M200 and c200.
    #the subroutine corenfw_mass will write out rs and rhos while
    #doing this:
    ranal = np.logspace(-1,3,500)
    manal = corenfw_mass(ranal,M200,c200,\
                         oden,tSF,Rhalf,rhocrit,eta,kappa)

    #Make a plot of the above profile:
    figsize = 8
    figx = figsize
    figy = figsize
    myfontsize = 30
    mylinewidth = 3
    mylinewidth2 = 4

    #Use LaTeX Fonts: 
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')

    #Set thick axes: 
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)

    #Make nice tick marks:
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    plt.loglog()

    plt.xlabel(r'$r\,[\mathrm{kpc}]$',fontsize=myfontsize)
    plt.ylabel(r'$M_\mathrm{DM}(<r)\,[\mathrm{M}_\odot]$',\
                   fontsize=myfontsize)

    plt.plot(ranal,manal,\
             linewidth=mylinewidth,color='black',\
             zorder=0)
                    
    plt.ylim([1e6,3e12])
    plt.xlim([0.1,300])
    
    plt.savefig('NFW.pdf')
