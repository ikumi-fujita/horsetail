########## HORSETAIL SIMULATION, R version ##########

# 'HTsim()' returns a matrix containing positions of the SPB-nucleus and MT tips.
# Execute 'Test.run()' to draw a graph of a simulation with the default parameters,
# which can be changed by giving a value to each argument: 'l', 'MTnumber',
# 'v' (as a vector), 'f' (as a vector), 'r', 'eta', 'Fpush', 'Fpull', 'dT', 'rs' and 'min'.

HTsim <- function( l=14.0,                # l, cell length [um]
                   MTnumber=3,            # MTnumber, number of MT on each side of the SPB
                   v=c(3.3,3.3,4.2,3.7),  # v, MT velocity [um/min]
                   # v[1],Vplus(cyt); v[2],Vplus(ctx); v[3],Vminus(cyt); v[4],Vminus(ctx)
                   f=c(0.001,0.03,0,0),   # f, catastroph or rescue frequency [/s]
                   # f[1],fcat(cyt); f[2],fcat(ctx); f[3],fres(cyt); f[4],fres(ctx)
                   r=2.0,                 # r, Stokes radius of the nucleus [um]
                   eta=1.0,               # eta, viscosity of cytosol [kg/m s]
                   Fpush=4.0,             # Fpush [pN]
                   Fpull=2.0,             # Fpull [pN]
                   dT=1,                  # dT, one step of simultion [s]
                   rs=60,                 # output results at every 'rs' steps
                   min=240                # min, duration of simulation [min]
) {
  
  ###### Parameters ######
  ### Parameters for simulation
  CellLength <- c(-l/2,l/2)               # CellLength, range of the cytoplasm
  InitialPosition = 0                     # IntialPosition, initial position of the SPB
  ST = 60*min/dT                          # ST, total steps
  
  ### Parameters for MT
  MTN2 <- MTnumber*2                      # MTN2, total number of MT
  MTNV <- c(1:MTN2)                       # an array of 1, 2,..., MTN2
  S <- rep(c(1,-1),length=MTN2)           # Direction of MT: 1,-1,1,-1,...
  Vplus <- c(v[1],v[2])/60*dT             # [um/step]
  Vminus <- c(v[3],v[4])/60*dT            # [um/step]
  Vmode <- c(Vplus[1],-Vminus[1],Vplus[2],-Vminus[2],0) # Vmode[MTmode].
                                          # When MT pushes the cortex (MTmode==5), MT length does not change.
  fcat <- c(f[1],f[2])*dT                 # fcat [/step]
  fres <- c(f[3],f[4])*dT                 # fres [/step]
  
  ### Parameters for forces
  Cdrag = 6.0*pi*r*eta                    # drag constant [kg/s]
  Cbuckle = 25*pi^2                       # buckling constant [pN um^2]
  PushLength = (Cbuckle/Fpush)^(1/2)      # If the length of an MT is longer than this length, the MT buckles.
  
  ###### Variables ######
  Position <- InitialPosition             # Position of the SPB
  MTlength <- numeric(MTN2)               # Length of MTs. Initialized to 0.
  MTposition <- rep(Position,MTN2)        # Position of MT tips. Initialized to the initial SPB position.
  MTmode <- rep(1,MTN2)                   # Modes of MTs. Initialized to 1.
  # 1: Growing in the cytoplasm
  # 2: Shrinking in the cytoplasm
  # 3: Growing along the cortex
  # 4: Shrinking along the cortex and generating pulling force
  # 5: Pushing the cortex
  
  ###### Initialization ######
  # Prepare a matrix to output the results (position of each MT at every 'rs' step)
  Result <- matrix(c(seq(0,ST/60,dT/60*rs),numeric((ST/rs+1)*(1+MTN2))),
                   ncol=(2+MTN2),nrow=(ST/rs+1),byrow=FALSE)
  colnames(Result) <- c("t","SPB",paste("MT",MTNV,sep=""))
  Result[1,2] = Position
  Result[1,3:(2+MTN2)] = MTposition
  
  ###### Run Simulation ######
  for (t in 1:ST) {
    ### Force Generation
    Force = 0
    for (x in MTNV) {
      # Pulling force
      if (MTmode[x]==4) {
        Force = Force + S[x]*Fpull
      }
      # Pushing force
      if (MTmode[x]==5) {
        Force = Force - S[x]*Fpush
      }
    }
    
    ### SPB movement
    # Change the length of each MT according to the MTmode
    MTlength <- MTlength + Vmode[MTmode]
    # Change the SPB position according to the Force
    Position = Position + Force/Cdrag*dT
    if (Position < CellLength[1]) {
      Position = CellLength[1]
    } else if (Position > CellLength[2]) {
      Position = CellLength[2]
    }
    # Change the tip position of each MT
    MTposition <- Position + S*MTlength
    
    for (x in MTNV) {
      ### Mode change depending on the position
      if (MTmode[x]==1) { # Growing in the cytoplasm
        if ((MTposition[x] <= CellLength[1])||(MTposition[x] >= CellLength[2])) {
          if (MTlength[x] < PushLength) {
            MTmode[x] = 5 # changed to Pushing the cortex
          } else {
            MTmode[x] = 3 # changed to Growing along the cortex
          }
        }
      }
      else if (MTmode[x]==2) { # Shrinking in the cytoplasm
        if ((MTposition[x] <= CellLength[1])||(MTposition[x] >= CellLength[2])) {
          MTmode[x] = 4 # changed to Shrinking along the cortex
        }
      }
      else if (MTmode[x]==3) { # Growing along the cortex
        if (MTlength[x] < PushLength) {
          MTmode[x] = 5 # changed to Pushing the cortex
        }
        if ((MTposition[x] > CellLength[1])&&(MTposition[x] < CellLength[2])) {
          MTmode[x] = 1 # changed to Growing in the cytoplasm
        }
      }
      else if (MTmode[x]==4) { # Shrinking along the cortex
        if ((MTposition[x] > CellLength[1])&&(MTposition[x] < CellLength[2])) {
          MTmode[x] = 2 # changed to Shrinking in the cytoplasm
        }
      }
      else if (MTmode[x]==5) { # Pushing the cortex
        if ((MTposition[x] > CellLength[1])&&(MTposition[x] < CellLength[2])) {
          MTmode[x] = 1 # changed to Growing in the cytoplasm
        }
      }
      
      ### catasrophe or rescue of each MT
      if (MTmode[x]==1) { # Growing in the cytoplasm
        if (runif(1)<fcat[1]) {
          MTmode[x] = 2 # changed to Shrinking in the cytoplasm
        }
      }
      else if (MTmode[x]==2) { # Shrinking in the cytoplasm
        if ((MTlength[x]<0.1)||(runif(1)<fres[1])) {
          MTmode[x] = 1 # changed to Growing in the cytoplasm
        }
      }
      else if (MTmode[x]==3) { # Growing along the cortex
        if (runif(1)<fcat[2]) {
          MTmode[x] = 4 # changed to Shrinking along the cortex
        }
      }
      else if (MTmode[x]==4) { # Shrinking along the cortex
        if ((MTlength[x]<0.1)||(runif(1)<fres[2])) {
          MTmode[x] = 3 # changed to Growing along the cortex
        }
      }
      else if (MTmode[x]==5) { # Pushing the cortex
        if (runif(1)<fcat[2]) {
          MTmode[x] = 4 # changed to Shrinking along the cortex
        }
      }
    }
    
    ### Output results once in every 'rs' step
    if ((t%%rs)==0) {
      Result[t/rs+1,2] <- Position
      Result[t/rs+1,MTNV+2] <- MTposition
    }
  }
  return(Result)
}

###### Run the simulation and draw a graph ######
Test.run <- function( l=14.0,
                      MTnumber=3,
                      v=c(3.3,3.3,4.2,3.7),
                      f=c(0.001,0.03,0,0),
                      r=2.0,
                      eta=1.0,
                      Fpush=4.0,
                      Fpull=2.0,
                      dT=1,
                      rs=60,
                      min=240
) {
  sim <- HTsim(l,MTnumber,v,f,r,eta,Fpush,Fpull,dT,rs,min)
  plot(sim[,1],sim[,2],type="l",ylim=c(-1*l,l),xaxs="i",yaxs="i",yaxp=c(-1*l,l,4),
       xlab=expression("Time (min)"),ylab=expression("x ("*mu*"m)"),lwd=3)
  for (i in 1:(ncol(sim)-2)) {
    lines(sim[,1],sim[,i+2],col=colors()[i*5+20],lwd=1)
  }
  abline(h=c(-l/2,l/2),lwd=0.5)
}
