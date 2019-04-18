#2018-2019 Joshua Fitzgerald
#Finds a trajectory matching a provided itinerary.
#Shapely, after some experimentation, has been abandoned as unreliable;
#manifold intersections are once again being found manually.
#Call this command as
#sudo python toboldlygo3.py energy mu itinerary
#for example (Jupiter) sudo python toboldlygo3.py -1.515 9.537e-4 x21
#If Python raises the exception RuntimeError: could not create GdkCursor object when the program
#is run or imported, you may need to run "export DISPLAY=:0.0" from the command line

##############################################################################################
#Utility Functions (code not called by the main program)                                     #
##############################################################################################    

def mufinder(m1, m2):
    """This function calculates mu given 2 masses. Don't worry about preventing m2 > m1;
       the function can account for that possibility."""

    if m2 > m1:
        m1, m2 = m2, m1

    return m2 / (m1 + m2)

def lowenergyinterval(mu):
    """This function finds the E2 and E3 values for given mu.
       Use it to choose an energy value."""

    xl2, _ = findx(mu, 2) #We find the x-coordinates of the Lagrange points 
    xl3, _ = findx(mu, 3)

    return findenergy(xl2, 0, 0, 0, mu), findenergy(xl3, 0, 0, 0, mu) #We then find the energy

    
##############################################################################################
#Core Trajectory Derivation Functions and Code                                               #
##############################################################################################    

def xroot(x, mu):
    """The equation of which we must find the root."""
    return -x + (mu * (-1 + mu + x))/abs(-1 + mu + x)**3 - ((-1 + mu)*(mu + x))/abs(mu + x)**3

def findx(mu, lnum):
    """Obtains the Hill sphere and x-coordinate for a mu-value and lnum."""
    
    hill = (mu/3)**(1.0/3.0)
    
    if lnum == 1: #lnum is used to request one of the collinear Lagrange points.
        guess = 1 - mu - hill * (1 - (1.0/3.0) * hill - (1.0/9.0) * hill ** 2)
    elif lnum == 2:
        guess = 1 - mu + hill * (1 + (1.0/3.0) * hill - (1.0/9.0) * hill ** 2)
    elif lnum == 3:
        guess = -1 #I know this isn't the formal guess the Mission Handbook might prescribe, but it should suffice
                   #as the L3 Lagrange point is the only collinear point with x < 0
    else:
        return "Invalid"

    return optimize.fsolve(xroot, guess, mu, xtol = 0.0)[0], hill #Solving so that the relative error between iterates seems to
                                                                 #be 0 should give machine precision

def findparams(mu, lnum, amplitude):
    """Finds the parameters required for stablizing a periodic orbit."""
    lagrangex, hill = findx(mu, lnum)
    if lagrangex == "Invalid":
         print "This application is only designed for use with the L1 or L2 points."
         sys.exit()
    
    barmu = mu * abs((lagrangex - 1.0 + mu)**-3.0) + (1.0 - mu) * abs((lagrangex + mu)**-3.0)
    nu = np.sqrt(-0.5 * (barmu - 2 - np.sqrt(9 * barmu ** 2.0 - 8 * barmu)))
    tau = -1.0 * (nu ** 2.0 + 2 * barmu + 1.0) / (2.0 * nu)
    vynew = -1.0 * amplitude * nu * tau

    print "----Calculated parameters----"
    print "Hill radius:", hill
    print "L" + str(lnum), "x-coordinate:", lagrangex
    print "Barred mu:", barmu, "\nNu:", nu, "\nTau:", tau, "\nStarting y-velocity:", vynew
    return lagrangex, vynew

def findenergy(x, y, vx, vy, mu):
    """Finds the energy at a point in the 4D phase-space."""
    r1 = ((x + mu)**2 + y**2) ** 0.5 
    r2 = ((x - (1 - mu))**2 + y**2) ** 0.5
    return 0.5 * (vx**2 + vy**2 - x**2 - y**2) - (1 - mu)/r1 - mu/r2 - 0.5 * (1 - mu) * mu

def uxx(x, y, m1, m2):
    return -1 - (3*m2*(-m1 + x)**2)/((-m1 + x)**2 + y**2)**(5.0/2.0) \
    + m2/((-m1 + x)**2 + y**2)**(3.0/2.0) - (3*m1*(m2 + x)**2)/((m2 + x)**2 + y**2)**(5.0/2.0) \
    + m1/((m2 + x)**2 + y**2)**(3.0/2.0)

def uxy(x, y, m1, m2):
    return -((3*m2*(-m1 + x)*y)/((-m1 + x)**2 + y**2)**(5.0/2.0)) \
           - (3*m1*(m2 + x)*y)/((m2 + x)**2 + y**2)**(5.0/2.0)

def uyy(x, y, m1, m2):
    return -1 - (3*m2*y**2)/((-m1 + x)**2 + y**2)**(5.0/2.0) \
    + m2/((-m1 + x)**2 + y**2)**(3.0/2.0) \
    - (3*m1*y**2)/((m2 + x)**2 + y**2)**(5.0/2.0) + m1/((m2 + x)**2 + y**2)**(3.0/2.0)

def rotating(x, y, t):
    return x * np.cos(t) + y * np.sin(t), -x * np.sin(t) + y * np.cos(t)

def rotatingv(x, y, vx, vy, t):
    return vx * np.cos(t) + vy * np.sin(t) + y, -vx * np.sin(t) + vy * np.cos(t) - x
    #return (vx - y) * np.cos(t) + (vy + x) * np.sin(t), (y - vx) * np.sin(t) + (vy + x) * np.cos(t)

def jacobian(t, refstate, m1, m2):
    """Finds the jacobian matrix at time t for the reference trajectory trajectory."""

    coordinates = refstate.transpose() #Load the trajectory at t

    uxx1 = uxx(refstate[0], refstate[1], m1, m2)
    uxy1 = uxy(refstate[0], refstate[1], m1, m2)
    uyy1 = uyy(refstate[0], refstate[1], m1, m2)
    return np.matrix([[0,     0,      1, 0],
                        [0,     0,      0, 1],
                        [-uxx1, -uxy1,  0, 2],
                        [-uxy1, -uyy1, -2, 0]])

def euler(t, state, refstate, h, m1, m2):
    """Use euler's method to integrate the state transition matrix."""
    return state + h * jacobian(t, refstate, m1, m2) * state

def RK4(t, state, refstate, h, m1, m2):
    """Use the fourth-order Runge-Kutta method to integrate the state transition matrix."""
    j = jacobian(t, refstate, m1, m2)
    F1 = h * j * state
    F2 = h * j * (state + 0.5 * F1)
    F3 = h * j * (state + 0.5 * F2)
    F4 = h * j * (state + F3)
    return state + (1.0/6.0) * (F1 + 2 * F2 + 2 * F3 + F4)

def simcycle(sim, positions, time):
    """Does the position-velocity integration at the next timestep."""
    sim.integrate(time) #Do the integration
    x1, y1 = rotating(sim.particles[2].x, sim.particles[2].y, time) #Get the rotating frame trajectory positions and velocities
    vx1, vy1 = rotatingv(x1, y1, sim.particles[2].vxyz[0], sim.particles[2].vxyz[1], time)
    positions[time] = np.array([x1, y1, vx1, vy1]) #Store them
    return sim, positions

def integratecycle(sim, mu, positions, mintegration, lasttime, timestep, solver, onlycurrentpos = False):
    """Calculates a trajectory at the next timestep."""
    time = lasttime + timestep
    sim, positions = simcycle(sim, positions, time)
    if mintegration is not None:
        mintegration[time] = solver(lasttime, mintegration[lasttime], positions[lasttime],
                                      timestep, 1 - mu, mu) #Integrate the state transition matrix
        
    if onlycurrentpos: #onlycurrentpos discards the positions it received and replaces them with what was calculated
        return sim, {time : positions[time]}, mintegration, time
    return sim, positions, mintegration, time

def integratecore(mu, timestep, x, y, vx, vy, conditionfun, args, findmatrix = True, starttime = 0.0):
    """Calculates a trajectory and state transition matrix until the conditionfun function is true. args go to conditionfun.
Coordinates must be in rotating frame. REBOUND uses the inertial frame, but this function outputs in the rotating frame (via
integratecycle)."""

    varlist = {"timestep":timestep} #We do this so that the crossing function can modify the timestep and have other saved-state variables
    
    sim = rebound.Simulation()
    sim.dt = 1e-14
    
    sim.add(m = 1 - mu) #Denormalized coordinates!
    sim.add(m = mu, a = 1.0)
    sim.move_to_com()
    sim.add(m = 0.0, x = x, y = y, vx = vx - y, vy = vy + x) #Does rotating to inertial conversion on velocity. Position coincides
                                                             #At t = 0 (which we assume)
    sim.move_to_com()
    
    ps = sim.particles       # ps is now an array of pointers and will change as the simulation runs

    #positions structure:
    #{timesa1:np.array(xa1, ya1, dxa1, dya1), timesa2:np.array(xa2, ya2, dxa2, dya2), ...}
    #positions is a dictionary for the particle. The keys are the elements of times; they access corresponding lists of position-velocity values.

    lasttime = starttime

    #mintegration: The matrix integration
    if findmatrix:
        mintegration = {0.0 : np.identity(4)}
    else:
        mintegration = None

    x1, y1 = rotating(ps[2].x, ps[2].y, lasttime) #Get the rotating frame trajectory positions and velocities
    vx1, vy1 = rotatingv(ps[2].x, ps[2].y, ps[2].vxyz[0], ps[2].vxyz[1], lasttime)
    positions = integratecycle(sim, mu, {0.0 : np.array([x1, y1, vx1, vy1])},
                               mintegration, 0.0, lasttime, euler, onlycurrentpos = True)[1] #We integrate to lasttime

    while True:
        o = integratecycle(sim, mu, positions, mintegration, lasttime, varlist["timestep"], euler)
        sim = o[0]
        positions = o[1]
        mintegration = o[2]
        time = o[3]

        lasttime += varlist["timestep"]
        
        #Conditionfun must first return a truth value. Conditionfun also specifies what integratecore will output, so
        #you must also return at least None. Remember to account for whether you use mintegration;
        #if you use mintegration in conditionfun, don't specify findmatrix = False, or it will constantly be
        #None.
        #Take advantage of the fact that conditionfun can edit varlist; you can change the timestep or save
        #state between tests using varlist.
        
        condout = conditionfun(mu, sim, positions, mintegration, time, varlist, args) 
        if condout[0]:
            break
        
    return condout[1:]
    
    
def crossingtest(mu, sim, positions, mintegration, time, varlist, vynew):
    """The test function to see whether we've hit the crossing in differential correction."""
    vynew = vynew[0] #Just in case multiple arguments are passed, we get the first
    if abs(positions[time][1]) < 1e-11: #If the last particle y value was at the crossing, break
        print "crossing point at time:", time
        print "x-coordinate:", positions[time][0]
        print "y-coordinate:", positions[time][1]
        print "x-velocity:", positions[time][2]
        print "y-velocity:", positions[time][3]
        print "State transition matrix at most recent t:\n", mintegration[time]
        vynew = vynew - positions[time][2] / (mintegration[time][2, 3] -
                positions[time][2] * mintegration[time][1, 3] / positions[time][3])
        print "Corrected y-velocity:", vynew
        return True, vynew, positions[time][2], time * 2.0 #We return the corrected vy, vx1, period values
    
    elif positions[time][1] * np.sign(varlist["timestep"]) < 0.0: #If we have changed sign but not entered the crossing zone, we alter
                                                          #the timestep
        varlist["timestep"] /= -2.0
        return False, None
    return False, None

def periodtest(mu, sim, positions, mintegration, time, varlist, period):
    """The test function to see whether the time is greater than or equal to the periodic orbit's period."""
    period = period[0]
    if time <= varlist["timestep"]:
        print "Start:", positions[0.0]
        print "Time: 0.0"
    elif time % (period / 10.0) < varlist["timestep"]:
        print positions[time]
        print "Time:", time
    elif time == period:
        print "End:", positions[time]
        print "Time:", time
        return True, mintegration[time]
    elif (time - period) * np.sign(varlist["timestep"]) > 0:
        varlist["timestep"] /= -2.0
        return False, None
    return False, None

def u1(mu, positions, time, timestep, energy, lastpos): return u14(mu, positions, time, timestep, energy, lastpos, 1) #Some wrappers that make things slightly easier
def u2(mu, positions, time, timestep, energy, lastpos): return u23(mu, positions, time, timestep, energy, lastpos, 2)
def u3(mu, positions, time, timestep, energy, lastpos): return u23(mu, positions, time, timestep, energy, lastpos, 3)
def u4(mu, positions, time, timestep, energy, lastpos): return u14(mu, positions, time, timestep, energy, lastpos, 4)

def venergy(x, y, v, energy, mu):
    """Interestingly, the functions for vx and vy obtained from the energy equation are nearly identical,
so we can just supply the one we've got and get the other."""
    r1 = np.sqrt((x + mu)**2 + y**2) 
    r2 = np.sqrt((x - (1 - mu))**2 + y**2)
    baru = -0.5 * (x**2 + y**2) - (1 - mu) / r1 - mu / r2 - 0.5 * (1 - mu) * mu
    return np.sqrt(-v ** 2 - 2 * baru + 2 * energy) #We return the positive value; we can change it later.

def u14(mu, positions, time, timestep, energy, lastpos, number):
    """The test for the U1 and U4 Poincare sections. Number is which section to use; specify 1 or 4."""

    cposition = positions[time]
    
    if number == 1:
        if cposition[0] <= -1.0: #We're trying to make sure that 1 doesn't get hits in the exterior realm
            return 0
        xlessthan = 0.0
    elif number == 4:
        xlessthan = -1.0
    else:
        raise ValueError
    if cposition[0] < xlessthan: #No point in doing anything else if we don't meet basic prerequisites.
        energyvy = venergy(cposition[0], cposition[1], cposition[2], energy, mu)
        #To be honest, I'm not really sure why the Mission Handbook has us finding vy from the
        #energy equation (I've got vy already stored as cposition[3]) but I'll get it
        #this way for now.
        
        if energyvy > 0.0:    
            if abs(cposition[1]) < 1e-11: #If y is close enough to 0
                return 2
            if np.sign(cposition[1]) != np.sign(lastpos[1]): #In other words, if we just crossed the section but are not close enough to 0
                return 1
    return 0

def u23(mu, positions, time, timestep, energy, lastpos, number):
    """The test for the U2 and U3 Poincare sections. Number is which section to use; specify 2 or 3."""
    cposition = positions[time]
    
    if number == 2:
        ydirection = -1.0
    elif number == 3:
        ydirection = 1.0
    else:
        raise ValueError

    if cposition[1] * ydirection > 0.0: #No point in doing anything else if we don't meet basic prerequisites.
        energyvx = venergy(cposition[0], cposition[1], cposition[3], energy, mu)
        #Again, unsure about this.
        if energyvx > 0.0:
            
            if abs(cposition[0] - (1.0 - mu)) < 1e-11: #If x is close enough to 1 - mu
                return 2
            if np.sign(cposition[0] - (1.0 - mu)) != np.sign(lastpos[0] - (1.0 - mu)): #In other words, if we just crossed the section
                                                                   #but are not close enough to 1 - mu
                return 1
    return 0

def manifoldtest(mu, sim, positions, mintegration, time, varlist, args):
    """The test function to get the manifolds with a Poincare section."""
    #args = [section, energy, maxexectime, depth, tolerance]
    #section is a function that takes a position-velocity point in the phase space and
    #checks to see whether we've crossed or are at the section.
    #maxexectime is the number of seconds to spend trying to find the trajectory before skipping it and moving on.
    #depth - 1 is how many times to ignore crossings before attempting a hit. This parameter allows us to find trajectories that wind around several times.
    #tolerance is how close the energy at the crossing must be to the trajectory energy 
    #If we have crossed the section but are not yet within 1e-11 of it,
    #we perform a correction to the timestep to get us closer.
    #If section returns a 0, we need to keep integrating. If it returns a 1, we need to correct the timestep.
    #If it returns a 2, we are sufficiently close to the section.
    #We need the energy for the section tests.

    if "lastpos" not in varlist:
        varlist["lastpos"] = positions[sorted(positions.keys())[-2]]

    if "starttime" not in varlist:
        varlist["starttime"] = systemtime.clock()

    if "remainhits" not in varlist:
        varlist["remainhits"] = args[3]

    if systemtime.clock() - varlist["starttime"] > args[2]:
        print "skipped"
        return True, None
    
    check = args[0](mu, positions, time, varlist["timestep"], args[1], varlist["lastpos"])

    varlist["lastpos"] = positions[time]

    if check and varlist["remainhits"] > 0: #If we have a crossing or a hit on the Poincare surface and remainhits is positive, we decrement remainhits
        varlist["remainhits"] -= 1
        print varlist["remainhits"]
    
    if check == 2 and varlist["remainhits"] <= 0: #If we have a hit and remainhits is not positive, we try to yield the crossing
        print "crossed"
        computedenergy = findenergy(positions[time][0], positions[time][1], positions[time][2], positions[time][3], mu)
        if abs(computedenergy - energy) < args[4]:
            #We ensure that the energy is within a specific tolerance (per Dr. Ross's suggestion); otherwise, we throw out the point
            return True, positions[time]
        print "Point energy", computedenergy, "not within tolerance interval for energy", energy, ". Point will be discarded"
        return True, None
        
    if check == 1 and varlist["remainhits"] <= 0: #If we have a crossing and remainhits is not positive, we change the timestep
            varlist["timestep"] /= -2.0
    
    return False, None

def visualizertest(mu, sim, positions, mintegration, time, varlist, args):
    """The test function for the visualizer."""
    #args = [endtime]; endtime is the time at which to stop integrating.
    if time >= args[0]:
        return True, positions, sim
    return False, None

def pocorrector(mu, lagrangex, timestep, amplitude, vynew):
    """Finds a correct Lagrange point periodic orbit for a mu-value, a Lagrange point x-coordinate, a timestep, an amplitude,
and a guess y velocity."""
    
    vx1req = False
    iterationcount = 1

    while True:
        print "----Beginning stabilization iteration " + str(iterationcount) + "----"

        vynew, vx1, period = integratecore(mu, timestep, lagrangex - amplitude, 0, 0, vynew, crossingtest, [vynew])
           
        if abs(vx1) < 1e-8:
            print "----Orbit stabilized for acceptable vx1. Use parameters from iteration " + str(iterationcount) + "----"
            break
        
        iterationcount += 1    

    return lagrangex - amplitude, vynew, period #We return the x-coordinate, y-velocity, and period

def energypo(energy, mu, lnum, timestep, tolerance, sa1, sa2):
    """Finds the peroidic orbit for an arbitrary energy and seed amplitudes sa1 and sa2."""
    print "FINDING SA1"
    lagrangex, sa1vynew = findparams(mu, lnum, sa1) #Find the seed orbits
    sa1c = pocorrector(mu, lagrangex, timestep, sa1, sa1vynew)
    print "FINDING SA2"
    _, sa2vynew = findparams(mu, lnum, sa2)
    sa2c = pocorrector(mu, lagrangex, timestep, sa2, sa2vynew)
    
    energies = [[findenergy(sa1c[0], 0, 0, sa1c[1], mu),sa1c], #Structure: [[energy1, [x1, vy1]], [energy2, [x2, vy2]], ...]
                [findenergy(sa2c[0], 0, 0, sa2c[1], mu),sa2c]] #We don't use a dict because we need everything sorted
    dindex = -2 #The index in energies used to calculate the delta values
    print "USING NUMERICAL CONTINUATION"
    while True:
        x = 2 * energies[-1][1][0] - energies[dindex][1][0]
        guessvy = 2 * energies[-1][1][1] - energies[dindex][1][1]
        x, vy, period = pocorrector(mu, lagrangex, timestep, lagrangex - x, guessvy)
        currentenergy = findenergy(x, 0, 0, vy, mu)
        print x, vy, period, currentenergy
        energies.append([currentenergy, [x, vy, period]])
        if currentenergy > energy: #If the current energy is greater than the desired energy
            if abs(currentenergy - energy) > tolerance:
                raise RuntimeError("Energy of periodic orbit not within specified tolerance")
            break
    print energies
    
    return x, vy, period

def monodromy(mu, x, vy, timestep, period):
    """Calculates the monodromy matrix for the periodic orbit with the initial condition
corresponding to x and vy with period period."""
    print "CALCULATING MONODROMY MATRIX"
    return integratecore(mu, timestep, x, 0, 0, vy, periodtest, [period])

def monodromyeigen(matrix):
    eigvals = np.linalg.eig(matrix)[0].real[0] #We get the eigenvectors that correspond to the correct eigenvalues
    maxeigval = max(eigvals)
    mineigval = min(eigvals)

    #1e-4: 1134.46960566 1137.0436898
    print "Eigenvalue test:", maxeigval, 1.0 / mineigval
    
    for count, i in enumerate(eigvals):
        if i == maxeigval:
            unstable = np.linalg.eig(matrix)[1].real[0][count]
        elif i == mineigval:
            stable = np.linalg.eig(matrix)[1].real[0][count]
            
    #raise RuntimeError("Breakpoint")

    return unstable, stable

def manifoldtype(itinerarypart, stability):
    """Finds the type of manifold needed."""
    #Finding the Lagrange point. itinerarypart will always contain a 2
    if "1" in itinerarypart:
        lnum = 1
    elif "x" in itinerarypart:
        lnum = 2

    if stability: #We use modifier to account for whether a manifold is stable or unstable
        modifier = 1.0
    else:
        modifier = -1.0

    #Finding the sign of the manifold branch.    
    if itinerarypart[0] == "1":
        sign = -modifier
    elif itinerarypart[0] == "x":
        sign = modifier
    elif itinerarypart[0] == "2":
        if itinerarypart[1] == "x":
            sign = -modifier
        elif itinerarypart[1] == "1":
            sign = modifier
            
    print lnum, sign
    return lnum, sign

def manifoldinitial(energy, mu, itinerarypart, timestep, mtimestep, tolerance,
                    seedorbit1, seedorbit2, stability, lowerdisplacement, upperdisplacement,
                    numdisplacements):
    """Finds the initial conditions for a manifold. itinerarypart
is a two-character string containing the corresponding two destinations."""

    lnum, sign = manifoldtype(itinerarypart, stability)

    x, vy, period = energypo(energy, mu, lnum, mtimestep, tolerance, seedorbit1, seedorbit2)
    
    mono = monodromy(mu, x, vy, mtimestep, period)
    unstable, stable = monodromyeigen(mono)

    print lowerdisplacement, upperdisplacement, numdisplacements

    epsilonvals = np.geomspace(sign * lowerdisplacement, sign * upperdisplacement, numdisplacements)
    manifoldregion = []
    
    if stability:
        eig = stable
    else:
        eig = unstable
        
    for i in epsilonvals:
        manifoldregion.append([x + i * eig[0], i * eig[1], i * eig[2], vy + i * eig[3]])

    return manifoldregion

def integrateregion(energy, mu, timestep, timeout, region, depth, tolerance, section):
    """Integrates a set of points denoted by region to section and returns the set of points at the depth-th intersection."""
    integratedregion = []
    for i in region:
        intpos = integratecore(mu, timestep, i[0], i[1], i[2], i[3], manifoldtest,
                                [section, energy, timeout, depth, tolerance], findmatrix = False)[0]
        if intpos is not None:
            integratedregion.append(intpos)
    return integratedregion

def createpolygon(region, section):
    """Creates and processes a manifold intersection list from a set of points."""
    
    if section == u1 or section == u4:
        points2D = ((i[0], i[2]) for i in region) #We assume that there is only one part to the input region 
    elif section == u2 or section == u3:
        points2D = ((i[1], i[3]) for i in region)

    return points2D

def getoverlap(energy, mu, region1, region2, section):
    """Finds the overlap between two regions provided as lists that contain coordinate pairs."""

    def plot_point(event, marker):
        if event.button == 2: #if only right clicks are valid for plotting points, the left and right mouse buttons 
            x = np.append(marker.get_xdata(), event.xdata) #are free to help with other things like zooming
            y = np.append(marker.get_ydata(), event.ydata)
            marker.set_xdata(x)
            marker.set_ydata(y)
            print event.xdata, event.ydata
            marker.figure.canvas.draw()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()

    print "-------OVERLAP VISUALIZER-------"
    print ("Visualizing tube overlap. Select points within the intersection region with the middle mouse button, " +
    "then close the interactive plot to continue. If you can't find any points in the intersection region, " +
    "close the plot without selecting anything and wait for more points to be generated. " +
    "Note that the first point selected will be the 'representative point' that will determine what trajectory is yielded if the " +
    "trajectory is complete; be sure that it is inside and not on the boundary of the intersection region.")

    #We plot the manifold intersection points and ask the user to find the tube overlap
    tb1x = []
    tb1y = []
    for i in region1:
        tb1x.append(i[0])
        tb1y.append(i[1])

    tb2x = []
    tb2y = []
    for i in region2:
        tb2x.append(i[0])
        tb2y.append(i[1])
    
    ax.plot(tb1x, tb1y, ".")
    ax.plot(tb2x, tb2y, ".")

    selected, = ax.plot([], [], ".")

    fig.canvas.mpl_connect("button_press_event",
                           lambda event : plot_point(event, selected))

    plt.show() #plt.show() is blocking; it waits for the plot window to be closed
    plt.close() #we do this even though the window was closed to be sure that everything was reset
     
    overlap = []

    xdata = selected.get_xdata() #It's worth noting that this "x" data may actually be y data
                                #because of the way the Poincare section coordinates are chosen
    ydata = selected.get_ydata() #This data will actually be for vx or vy

    if section == u1 or section == u3: #This logic allows us to have the correct
        vdirection = -1.0 #vx or vy value after solving for it as venergy returns positive values
    elif section == u2 or section == u4:
        vdirection = 1.0

    for i, _ in enumerate(xdata):
        
        if section == u1 or section == u4: #If we received a valid response, we convert the 2D coordinates back to the 4D phase space
            vy = vdirection * venergy(xdata[i], 0, ydata[i], energy, mu)
            if np.isfinite(vy):
                overlap.append([xdata[i], 0, ydata[i], vy]) #and store them
        elif section == u2 or section == u3:
            vx = vdirection * venergy(1.0 - mu, xdata[i], ydata[i], energy, mu)                         
            if np.isfinite(vx):
                overlap.append([1.0 - mu, xdata[i], vx, ydata[i]])
    
    print overlap
    
    if not overlap:

        print ("No overlap points selected. Retrying intersection:")
            
        return overlap, None, True

    return overlap, overlap[0], False
        

def itineraryisvalid(itinerary):
    """Verifies that an itinerary is able to be processed."""
    itinerary = itinerary.lower()
    if len(itinerary) < 3: #We don't allow itineraries that don't have at least three destinations for now
        return False
    for count, i in enumerate(itinerary):
        if i not in ("1","2","x"): #We only allow these three symbols
            return False
        if count + 1 < len(itinerary): #We don't allow travel between the m1 and exterior realms
                                       #or the same destination twice
            if i == "1" and itinerary[count + 1] == "x":
                return False
            if i == "x" and itinerary[count + 1] == "1":
                return False
            if i == itinerary[count + 1]:
                return False
    return True

def findsection(founditin, newitin):
    """Finds the appropriate Poincare section for given itinerary parts."""
    if founditin[-1] == "1":
        return u1
    if founditin[-1] == "x":
        return u4
    if founditin[-1] == "2":
        if newitin[1] == "1":
            return u3
        if newitin[1] == "x":
            return u2
    
def posreal(value, exception = ValueError):
    """Casts to positive real (well, float) if possible"""
    newvalue = float(value) #The exception may get raised here by the float function if needed
    if newvalue <= 0:
        raise exception
    return newvalue

def natural(value):
    """Casts to natural (integer) if possible"""
    newvalue = int(value) #The exception may get raised here by the int function if needed
    if newvalue <= 0:
        raise exception
    return newvalue

def whole(value):
    """Casts to whole (integer) if possible"""
    newvalue = int(value) #The exception may get raised here by the int function if needed
    if newvalue < 0:
        raise exception
    return newvalue

def posrealargparse(value): return posreal(value, argparse.ArgumentTypeError)
def naturalargparse(value): return natural(value, argparse.ArgumentTypeError)
def wholeargparse(value): return whole(value, argparse.ArgumentTypeError)
    
def getinterceptparams(initialinterceptparams, interceptparams):
    """Gets parameters for an upcoming interception attempt from the user."""
    #initialinterceptparams are default values of the parameters without the prefix
    #and for the corresponding manifolds

    #We return a new version of interceptparams with any requested changes made

    #The structure of interceptparams lists is
    #[stablowerdisplacement, stabupperdisplacement, stabnumdisplacements,
    #unstablowerdisplacement, unstabupperdisplacement, unstabnumdisplacements, depth]

    #We also return a list, recommendations, with two Boolean values.
    #Each value recommends whether we should recompute the stable or unstable, respectively, manifolds based on what values changed

    print "\n-------PARAMETER EDITING PROMPT-------"

    print "\nCURRENT VALUES:"
    print "\n\tFOR DISPLACEMENTS OF STABLE MANIFOLD INITIAL CONDITIONS:"
    print "\t\tslb - the lower bound: " + str(interceptparams[0])
    print "\t\tsub - the upper bound: " + str(interceptparams[1])
    print "\t\tsn - the number of displacements: " + str(interceptparams[2])
            
    print "\n\tFOR DISPLACEMENTS OF UNSTABLE MANIFOLD INITIAL CONDITIONS:"
    print "\t\tulb - the lower bound: " + str(interceptparams[3])
    print "\t\tuub - the upper bound: " + str(interceptparams[4])
    print "\t\tun - the number of displacements: " + str(interceptparams[5])

    print "\n\tFOR MANIFOLD INTERCEPTION DEPTH:"
    print "\t\tdepth - alter the desired depth of the manifold interception: " + str(interceptparams[6])
    
    print "\nPARAMETER EDITING PROMPT INSTRUCTIONS:"
    print "\tAll parameter change commands are of the form [command] [value]."
    print "\tYou can also type default to get a list of default parameters you can use."
    print "\tWhen you are done, type done"

    print "\nFOR DISPLACEMENTS OF STABLE MANIFOLD INITIAL CONDITIONS:"
    print "\tslb - alter the lower bound (positive real number required for value argument)"
    print "\tsub - alter the upper bound (positive real number required for value argument)"  
    print "\tsn - alter the number of displacements (natural number required for value argument)" 
    
    print "\nFOR DISPLACEMENTS OF UNSTABLE MANIFOLD INITIAL CONDITIONS:"
    print "\tulb - alter the lower bound (positive real number required for value argument)"
    print "\tuub - alter the upper bound (positive real number required for value argument)"  
    print "\tun - alter the number of displacements (natural number required for value argument)"

    print "\nFOR MANIFOLD INTERCEPTION DEPTH:"
    print "\tdepth - alter the desired depth of the manifold interception (natural number required for value argument)"

    commands = {"slb":[0,posreal,1],  #this dict maps between commands for parameters, corresponding parameter array placements, functions
                "sub":[1,posreal,1],  #that try to handle the value parameters inside commands (to cast them to acceptable types
                "sn":[2,natural,1],   #or to reject them), and affiliation (see below). "R+" means positive real, "J" means natural number, and "W" means whole number
                "ulb":[3,posreal,-1],
                "uub":[4,posreal,-1], #affiliation is what manifold a value is associated with. It is used in making recomputation recommendations.
                "un":[5,natural,-1],  #0 means no affiliation, 1 means stable manifold affiliation, and -1 means unstable manifold affiliation
                "depth":[6,natural,0]}
                                      #and yes, I know that at times like this a dedicated commandproperties object might make sense

    recommendations = [False, False] #no recommendations to recompute are made until changes to corresponding variables are made
        
    while True:
        resp = raw_input(">>> ").lower()
        if resp == "done":
            break
        if resp == "default":
            print "\nFOR DISPLACEMENTS OF STABLE MANIFOLD INITIAL CONDITIONS:"
            print "\tslb - the lower bound: " + str(initialinterceptparams[0])
            print "\tsub - the upper bound: " + str(initialinterceptparams[1])
            print "\tsn - the number of displacements: " + str(initialinterceptparams[2])
            
            print "\nFOR DISPLACEMENTS OF UNSTABLE MANIFOLD INITIAL CONDITIONS:"
            print "\tulb - the lower bound: " + str(initialinterceptparams[3])
            print "\tuub - the upper bound: " + str(initialinterceptparams[4])
            print "\tun - the number of displacements: " + str(initialinterceptparams[5])

            print "\nFOR MANIFOLD INTERCEPTION DEPTH:"
            print "\tdepth - alter the desired depth of the manifold interception: " + str(initialinterceptparams[6])
            continue
        resp = resp.split(" ")
        if len(resp) < 2:
            print "Please provide a valid response."
            continue
        if resp[0] in commands:
            try:
                #We get the value provided as the command argument and process it
                value = commands[resp[0]][1](resp[1])
            except ValueError:
                print "Argument type invalid."
                continue
            interceptparams[commands[resp[0]][0]] = value #We then store a key-value pair with the index from the command dict and the value
            if commands[resp[0]][2] == 1: #We also update the recommendation list
                recommendations[0] = True
                print "Note: the value of a parameter for the stable manifold has been modified. The initial condition region will be recomputed"
            if commands[resp[0]][2] == -1:
                recommendations[1] = True
                print ("Note: the value of a parameter for the unstable manifold has been modified. If the unstable manifold is still being integrated, " +
                "the initial condition region will be recomputed")
        else:
            print "Please provide a valid response."

    return interceptparams, recommendations

def finditinerary(energy, mu, itinerary, timestep, mtimestep, ttimestep, timeout, tolerance, seedorbit1, seedorbit2,
                  initialinterceptparams):
    """Finds an initial condition corresponding to an itinerary."""

    if not itineraryisvalid(itinerary):
        return "invalid itinerary"

    #founditin is the part of the itinerary for which we have a trajectory.
    founditin = itinerary[:2]

    findinginitial = True #findinginitial signifies that the initial unstable manifold can be found

    #forwardregion is the part of the trajectory to be integrated forward.
    forwardregion = []

    #The structure of the interceptparams list is as follows:
    #[unstablowerdisplacement, unstabupperdisplacement, unstabnumdisplacements,
    #stablowerdisplacement, stabupperdisplacement, stabnumdisplacements, depth].
    interceptparams = initialinterceptparams
    
    while founditin != itinerary:
        #We get the initial conditions for the stable manifold
        stableitin = itinerary[len(founditin)-1:len(founditin)+1] #stableitin is the itinerary for the stable manifold

        #We get the Poincare section we need
        section = findsection(founditin, stableitin)
        print section.__name__

        backwardregion = []

        unstabcomputemic = True #These variables determine whether the manifold initial conditions will be computed
        stabcomputemic = True

        while True: #This loop allows us to edit the properties of the upcoming integration

            forwardintegration = []
            backwardintegration = []

            interceptparams, recommendations = getinterceptparams(initialinterceptparams, interceptparams)
            print interceptparams

            if recommendations[0] or not backwardregion:
                print "-------FINDING STABLE MANIFOLD INITIAL CONDITIONS-------"
                backwardregion = manifoldinitial(energy, mu, stableitin, timestep, mtimestep, tolerance,
                                                 seedorbit1, seedorbit2, True, interceptparams[3], interceptparams[4],
                                                 interceptparams[5])

                stabcomputemic = False
            
            if (recommendations[1] and findinginitial) or not forwardregion: #If new manifold initial conditions are requested, we provide them
                print "-------FINDING UNSTABLE MANIFOLD INITIAL CONDITIONS-------"
                forwardregion = manifoldinitial(energy, mu, itinerary[:2], timestep, mtimestep, tolerance,
                                                seedorbit1, seedorbit2, False, interceptparams[0], interceptparams[1],
                                                interceptparams[2])
                unstabcomputemic = False
            
            print "-------Adding manifold intercepts at depth " + str(interceptparams[6]) + "-------"

            #We integrate the existing trajectory piece forward and the stable manifold backwards to the section
            print "-------INTEGRATING EXISTING TRAJECTORY FORWARDS-------"
            newforward = createpolygon(integrateregion(energy, mu, ttimestep, timeout, forwardregion, interceptparams[6], tolerance, section), section)
            forwardintegration += newforward #We get the new integration and add it to what we have
            print "-------INTEGRATING STABLE MANIFOLD BACKWARDS-------"
            newbackward = createpolygon(integrateregion(energy, mu, ttimestep * -1, timeout, backwardregion, interceptparams[6], tolerance, section), section)
            backwardintegration += newbackward
            
            #We find the new forward region as the overlap between the two regions
            forwardregiontemp, rep, stay = getoverlap(energy, mu, forwardintegration, backwardintegration, section)

            if not stay:
                
                forwardregion = forwardregiontemp
                break
        
        #depth += 1
            
        print forwardregion

        founditin = founditin + itinerary[len(founditin)]

        print "-------INTERCEPT REGION FOUND FOR ITINERARY SUBSET " + founditin + " AT " + str(interceptparams[6]) + "-------"

        findinginitial = False

    print "SAMPLE INITIAL CONDITION:", rep
    return rep

def trajectoryviewer(mu, timestep, initial):
    """A command-line generator for trajectories corresponding to initial conditions initial."""
    
    print "-------TRAJECTORY VISUALIZER-------"
    print "Please input visualizer options in the form '[starttime] [endtime]'. To quit, type 'done'. "
    #print "Note that the visualization will be saved as [filename].html."
    while True:
        options = raw_input(">>> ")
        options = options.strip()
        if options == "done":
            break
        else:
            options = options.split(" ")
            try:
                start = float(options[0]) #We try to get the values
                end = float(options[1])
                #filename = str(options[2])
            except:
                print "Please provide a valid response."
                continue

            print "Visualizing..."
            positions, sim = integratecore(mu, timestep, initial[0], initial[1], initial[2], initial[3],
                                           visualizertest, [end], findmatrix = False, starttime = start) 
            keys = sorted(positions.keys())

            fig = plt.figure() #Set up plot
            ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1.5, 1.5), ylim=(-1.5,1.5))
            ax.grid()
            ax.set_aspect("equal")

            viewpos = [] #The positions of particles in a form processable by the viewer

            for i in sim.particles:
                viewpos.append(ax.plot([], [], '.', markersize = 10.0)[0])
            

            def _init():
                """An init function for the visualizer."""
                for i in viewpos: 
                    i.set_data([], [])
                return viewpos

            def _animate(i):
                """An animation function for the visualizer."""
                if i == 0:
                    viewpos[2].set_data([], [])
                
                viewpos[0].set_data([-mu], [0]) #The positions of the masses stay constant
                viewpos[1].set_data([1.0 - mu], [0]) #So we can just keep giving the same value

                xlist = viewpos[2].get_data()[0]
                ylist = viewpos[2].get_data()[1]
                
                xlist.append(positions[keys[i]][0])
                ylist.append(positions[keys[i]][1])
                
                viewpos[2].set_data(xlist, ylist) #The position of the particle does not
                return viewpos

            ani = animation.FuncAnimation(fig, _animate,
                                          np.arange(0, len(keys)), interval=25, blit=True, init_func=_init)

            #ani.save(filename + ".html", dpi = 200, fps=60)

            plt.show()

def finditinvisual(energy, mu, itinerary, timestep, mtimestep, ttimestep, vtimestep, timeout, tolerance, seedorbit1,
                   seedorbit2, initialinterceptparams):
    """A version of finditinerary with the visualizer attached."""
    cond = finditinerary(energy, mu, itinerary, timestep, mtimestep, ttimestep, timeout, tolerance, seedorbit1,
                        seedorbit2, initialinterceptparams)
    if cond != "invalid itinerary":
        trajectoryviewer(mu, vtimestep, cond)
    else:
        print cond

import matplotlib
matplotlib.use("GTKAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as mpimg
import rebound
import sys
import time as systemtime
import numpy as np
from scipy import optimize
import random
import argparse

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "Finds a trajectory matching a provided itinerary")
    parser.add_argument("energy", type=float, help="the energy value of the desired trajectory")
    parser.add_argument("mu", type=float, help="the system's mass parameter")
    parser.add_argument("itinerary", type=str, help="the desired trajectory's itinerary")
    parser.add_argument("-t", "--timestep", type=float, help="the timestep for integrating anything without a dedicated timestep")
    parser.add_argument("-mt", "--monodromytimestep", type=float, help="the timestep for integrating monodromy matrices")
    parser.add_argument("-tt", "--tubetimestep", type=float, help="the timestep for integrating manifold/tube regions")
    parser.add_argument("-vt", "--visualizationtimestep", type=float, help="the timestep for integrating visualizations")
    parser.add_argument("-to", "--timeout", type=float, help="the maximum time to spend integrating a tube boundary")
    parser.add_argument("-tl", "--energytolerance", type=float, help="Energies must be this close to the prescribed energy")
    parser.add_argument("-so1", "--seedorbit1", type=float, help="a smaller seed orbit for numerical continuation")
    parser.add_argument("-so2", "--seedorbit2", type=float, help="a larger seed orbit for numerical continuation")
    parser.add_argument("-slb", "--stablelowerbound", type=posrealargparse, help="the default lower bound for the displacement for the stable manifold initial conditions")
    parser.add_argument("-sub", "--stableupperbound", type=posrealargparse, help="the default upper bound for the displacement for the stable manifold initial conditions")
    parser.add_argument("-sn", "--stablenumber", type=naturalargparse, help="the default number of displacements for the stable manifold initial conditions")
    parser.add_argument("-ulb", "--unstablelowerbound", type=posrealargparse, help="the default lower bound for the displacement for the unstable manifold initial conditions")
    parser.add_argument("-uub", "--unstableupperbound", type=posrealargparse, help="the default upper bound for the displacement for the unstable manifold initial conditions")
    parser.add_argument("-un", "--unstablenumber", type=naturalargparse, help="the default number of displacements for the unstable manifold initial conditions")
    parser.add_argument("-d", "--depth", type=naturalargparse, help="the default depth at which to integrate manifolds")

    inputstuff = parser.parse_args()

    energy = inputstuff.energy
    mu = inputstuff.mu
    itinerary = inputstuff.itinerary

    if inputstuff.timestep:
        timestep = inputstuff.timestep
    else:
        timestep = 1e-4

    if inputstuff.monodromytimestep:
        mtimestep = inputstuff.monodromytimestep
    else:
        mtimestep = 1e-3

    if inputstuff.tubetimestep:
        ttimestep = inputstuff.tubetimestep
    else:
        ttimestep = 1e-3

    if inputstuff.visualizationtimestep:
        vtimestep = inputstuff.visualizationtimestep
    else:
        vtimestep = 2e-2

    if inputstuff.timeout:
        timeout = inputstuff.timeout
    else:
        timeout = 10.0

    if inputstuff.energytolerance:
        tolerance = inputstuff.energytolerance
    else:
        tolerance = 1e-3

    if inputstuff.seedorbit1:
        seedorbit1 = inputstuff.seedorbit1
    else:
        seedorbit1 = 0.0001

    if inputstuff.seedorbit2:
        seedorbit2 = inputstuff.seedorbit2
    else:
        seedorbit2 = 0.0004

    initialinterceptparams = []

    if inputstuff.stablelowerbound:
        initialinterceptparams.append(inputstuff.stablelowerbound)
    else:
        initialinterceptparams.append(1e-8)

    if inputstuff.stableupperbound:
        initialinterceptparams.append(inputstuff.stableupperbound)
    else:
        initialinterceptparams.append(1e-4)

    if inputstuff.stablenumber:
        initialinterceptparams.append(inputstuff.stablenumber)
    else:
        initialinterceptparams.append(200)

    if inputstuff.unstablelowerbound:
        initialinterceptparams.append(inputstuff.unstablelowerbound)
    else:
        initialinterceptparams.append(1e-8)

    if inputstuff.unstableupperbound:
        initialinterceptparams.append(inputstuff.unstableupperbound)
    else:
        initialinterceptparams.append(1e-4)

    if inputstuff.unstablenumber:
        initialinterceptparams.append(inputstuff.unstablenumber)
    else:
        initialinterceptparams.append(200)

    if inputstuff.depth:
        initialinterceptparams.append(inputstuff.depth)
    else:
        initialinterceptparams.append(1) 

    finditinvisual(energy, mu, itinerary, timestep, mtimestep,
                    ttimestep, vtimestep, timeout, tolerance, seedorbit1, seedorbit2, initialinterceptparams)
