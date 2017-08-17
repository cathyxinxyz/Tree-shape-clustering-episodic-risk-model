
import random 
import numpy as np
import matplotlib.pyplot as plt
from random import Random, expovariate
import networkx as nx
import sys
import cPickle as pickle
from Parms import Parms
from Agent import Agent

MaxSimTime, simn, pickled_as=sys.argv[1:5]

def Initiate_Parms():    
    Parms.fL=1-Parms.fH
    Parms.B2=Parms.B1/Parms.r
    Parms.Ch=Parms.C*Parms.rho/(Parms.fL+Parms.rho*Parms.fH)
    Parms.Cl=Parms.C/(Parms.fL+Parms.rho*Parms.fH)
    Parms.Chh=Parms.Ch*Parms.mH
    Parms.Chg=Parms.Ch*(1-Parms.mH)
    Parms.Cll=Parms.Cl*Parms.mL
    Parms.Clg=Parms.Cl*(1-Parms.mL)

class Epid(object):
	#inherits from the generic object object
    def __init__(self):
        #define model paramters upon initilization         
        self.PPL={'Sh':[],'Sl':[],'I1h':[],'I1l':[],'I2h':[],'I2l':[]}
        self.Traj={'Sh':[],'Sl':[],'I1h':[],'I1l':[],'I2h':[],'I2l':[],'Time':[]}
        self.t = 0
        self.pplID=0
        self.data=[]         #record time, involved individuals and their states for events including transmission, disease progression, risk_change, or death  
        self.G=None          #complete transmission network that will be formulated as an networkx object
        self.Incidence=[]    #list of time and id of incident infections 
        self.Allevents=[]   
        self.Seed=None
        self.Sample_scenarios=None    
        self.Samples=list()
        self.mini_sample_t=None
        self.MaxSimTime=None
        
    def getID(self):
        myID=self.pplID
        self.pplID +=1
        return myID  

    def populate_model(self):
        #populate the model state with agents at t=0
        for state in ['Sh','Sl','I1h','I1l','I2h','I2l']:
            for i in range(int(Parms.initN[state])):                
                P=Agent(ID=self.getID(), state=state)
                P.remove_t=self.MaxSimTime
                self.PPL[state].append(P)   
                if 'I' in state:
                    self.Seed=P
                    P.transmit_time.append(0)
                    self.Incidence.append(P)
    def event_rates(self): # function that calculates the rates of each type of event that are used for stochastic simulation using Gillespie algorithm        
        #construct dictionary rates, with its keys as one of the following events:'birth','death','infection','change risk levels','disease progression'
        rates={}
        #birth events
        rates[('Sh', 'birth')] = Parms.N * Parms.m*Parms.fH
        rates[('Sl', 'birth')] = Parms.N * Parms.m*(1-Parms.fH)
        #death
        for state in ['Sh','Sl','I1h','I1l','I2h','I2l']:
            rates[(state,'death')] = len(self.PPL[state]) * Parms.m
	    #transmission events
        self.Lambdas(rates)
        #vaccination events
        for state in ['Sh','Sl','I1h','I1l', 'I2h', 'I2l']:
            if 'h' in state:
                rates[((state,state[:-1]+'l'), 'risk_change')] = len(self.PPL[state])*Parms.w*(1-Parms.fH)
            else:
                rates[((state,state[:-1]+'h'), 'risk_change')] = len(self.PPL[state])*Parms.w*Parms.fH
        
        #disease process
        for state in ['I1h','I1l', 'I2h', 'I2l']:
            if '1' in state:
                rates[((state,state.replace(state[1], '2')), 'state_change')]=len(self.PPL[state])*Parms.g1
            else:
                rates[((state,None), 'state_change')]=len(self.PPL[state])*Parms.g2

        return rates
        
    def Lambdas(self, rates): 
    ##function that calculates the rate of transmissions between infectives with speicific risk phase and infection state and susceptible individuals in a specific risk phase
    #denotations for infectivees: 'I1h':acutely infected in high risk phase, 'I2h':chronically infected in high risk phase, 'I1l':acutely infected in high risk phase, 'I2h':chronically infected in high risk phase, 'I2l': chronically infected in low risk phase:
    #denotations for susceptibles:'Sh': high-risk susceptible individuals, 'Sl': low-risk susceptible individuals
        Totcont_h=Parms.Ch*(len(self.PPL['Sh'])+len(self.PPL['I1h'])+len(self.PPL['I2h']))
        Totcont_g=Parms.Ch*(1-Parms.mH)*(len(self.PPL['Sh'])+len(self.PPL['I1h'])+len(self.PPL['I2h']))+Parms.Cl*(1-Parms.mL)*(len(self.PPL['Sl'])+len(self.PPL['I1l'])+len(self.PPL['I2l']))
        Totcont_l=Parms.Cl*(len(self.PPL['Sl'])+len(self.PPL['I1l'])+len(self.PPL['I2l']))
        fhh=len(self.PPL['Sh'])*Parms.Ch/Totcont_h
        fhg=len(self.PPL['Sh'])*Parms.Ch*(1-Parms.mH)/Totcont_g
        flg=len(self.PPL['Sl'])*Parms.Cl*(1-Parms.mL)/Totcont_g
        fll=len(self.PPL['Sl'])*Parms.Cl/Totcont_l
        Lambda_hh1=Parms.Ch*Parms.mH*len(self.PPL['I1h'])*Parms.B1            #
        Lambda_hh2=Parms.Ch*Parms.mH*len(self.PPL['I2h'])*Parms.B2
        Lambda_hg1=Parms.Ch*(1-Parms.mH)*len(self.PPL['I1h'])*Parms.B1
        Lambda_hg2=Parms.Ch*(1-Parms.mH)*len(self.PPL['I2h'])*Parms.B2
        Lambda_lg1=Parms.Cl*(1-Parms.mL)*len(self.PPL['I1l'])*Parms.B1
        Lambda_lg2=Parms.Cl*(1-Parms.mL)*len(self.PPL['I2l'])*Parms.B2
        Lambda_ll1=Parms.Cl*Parms.mL*len(self.PPL['I1l'])*Parms.B1
        Lambda_ll2=Parms.Cl*Parms.mL*len(self.PPL['I2l'])*Parms.B2
        rates[(('Sh', 'I1h'),'infection')]=Lambda_hh1*fhh+Lambda_hg1*fhg
        rates[(('Sh', 'I2h'),'infection')]=Lambda_hh2*fhh+Lambda_hg2*fhg
        rates[(('Sh', 'I1l'),'infection')]=Lambda_lg1*fhg
        rates[(('Sh', 'I2l'),'infection')]=Lambda_lg2*fhg
        rates[(('Sl', 'I1h'),'infection')]=Lambda_hg1*flg
        rates[(('Sl', 'I2h'),'infection')]=Lambda_hg2*flg
        rates[(('Sl', 'I1l'),'infection')]=Lambda_ll1*fll+Lambda_lg1*flg
        rates[(('Sl', 'I2l'),'infection')]=Lambda_ll2*fll+Lambda_lg2*flg


    def next_event(self, rates):        
        #determines the time and type of the next event given the event rates based on Gillespie algorithm
        total_rate = sum(rates.values()) 
        rnd_index = random.random() * total_rate #rnd float in [0, total_rate]
	
        #randomly draw time for next event from exponential distribution with mean as 1/total_rate
        self.t = self.t + np.random.exponential(scale=total_rate**-1)     
        
        #determine intervals of cumulative distribution of event rates
        rates_range={rates.keys()[n]:(sum(rates.values()[:n]), sum(rates.values()[:(n+1)])) for n in range(len(rates.values()))}
        self.Traj['Time'].append(self.t)
        
        #determine event types based on the interval of cumulative distribution that randomly drawn rate falls into
        for rr in rates_range.items():  
            if rnd_index>=rr[1][0] and rnd_index<rr[1][1]:
                event=rr[0]
        
                
        if 'birth' in event:
            state=event[0]
            ppl=Agent(ID=self.getID(), state=state)
            self.PPL[state].append(ppl)           
                     
        elif 'death' in event:
            state=event[0]
            ppl=self.PPL[state].pop(np.random.randint(0, len(self.PPL[state])))
            ppl.remove_t=self.t       #record death(removal) time
            self.data.append((self.t, ppl.ID, state, None, None, None))

        elif 'risk_change' in event:
            state=event[0][0]
            new_state=event[0][1]
            ppl = self.PPL[state].pop(np.random.randint(0, len(self.PPL[state])))
            ppl.state=new_state
            self.PPL[new_state].append(ppl)
            self.data.append((self.t, ppl.ID, state, new_state, None, None))
                      
        elif 'state_change' in event:
            state=event[0][0]
            new_state=event[0][1]
            ppl = self.PPL[state].pop(np.random.randint(0, len(self.PPL[state])))
            if new_state !=None:
                ppl.infectious_state=new_state
                self.PPL[new_state].append(ppl)               
            else:
                ppl.remove_t=self.t
            self.data.append((self.t, ppl.ID, state, new_state, None, None))
  
        elif 'infection' in event:         
            sus=self.PPL[event[0][0]].pop(np.random.randint(0, len(self.PPL[event[0][0]])))
            
            transmittor=random.sample(self.PPL[event[0][1]], 1)[0]
                       
            sus.state='I1'+sus.state[-1]
            self.PPL[sus.state].append(sus)
            self.data.append((self.t, sus.ID, event[0][0], sus.state, transmittor.ID, transmittor.state))

            #record the transmission history that is needed to construct transmission tree and graph
            sus.infectby=transmittor
            sus.state_at_infection=sus.state
            sus.infector_state=transmittor.state
            sus.infecttime=self.t
            transmittor.infectee.append(sus)
            sus.transmit_time.append(self.t)
            transmittor.infector=True
            
            sus.branchlength=self.t-transmittor.transmit_time[-1]
                       
            transmittor.transmit_time.append(self.t)            
                                                       
            self.Incidence.append(sus)
        
        #record the size of subpopualtion in each state at each time step
        for state in ['Sh','Sl','I1h','I1l','I2h','I2l']:
            self.Traj[state].append(len(self.PPL[state]))

           
#function that identify the instance by ID
def IdentifybyID(ID, epi):
    for P in epi.PplList:
        if P.ID==ID:
            return P 
            
def Transform_sto_traj(traj, scenario='hetero'):
    if scenario=='hetero':
        return np.array([[traj['Sh'][i],traj['Sl'][i],traj['I1h'][i],traj['I1l'][i],traj['I2h'][i],traj['I2l'][i]] for i in range(len(traj['Sh']))]) 
    elif scenario=='homo':
        return np.array([[traj['S'][i],traj['I1'][i],traj['I2'][i]] for i in range(len(traj['S']))])
        
        

def main(): 
    #MaxSimTime=10
    #sim_no=1
    #pickled_as='C:/git/transmission_tree_episodic_risk_model/epi_data.pkl'
    
    MaxSimTime=float(MaxSimTime)
    simn=float(simn)
    
    
    j=0
    while j<simn:           
        epi = Epid()  
        epi.MaxSimTime=MaxSimTime          
        epi.populate_model()
        
        while epi.t<MaxSimTime:
            rates = epi.event_rates()
            epi.next_event(rates)
        print epi.Incidence
        with open(pickled_as, 'wb') as output:
            pickle.dump(epi, output, pickle.HIGHEST_PROTOCOL)
            del epi
        j+=1
    
    
    
if __name__ == '__main__': main()
     
            
           
