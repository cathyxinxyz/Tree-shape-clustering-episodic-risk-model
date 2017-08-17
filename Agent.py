# -*- coding: utf-8 -*-

#A simple Agent class 
class Agent(object):	
    def __init__(self, ID, state):
        #init function not strictly nessisary here
        # given that the agent has no features
        self.event_id=0
        self.ID=ID
        self.state=state
        self.state_at_infection=None        #state at infection
        self.infector_state=None            #state of infector at time of infection
        self.infectby=None                  #the individual got disease from whom
        self.infectee=[]                    #list of infectees
        self.infector=False                 #whether the individual transmit disease to others      
        self.infecttime=None                #time of being infected
        self.transmit_time=[]         
        self.remove_t=float('inf')
        self.branchlength=None
        self.branch_other_events=None
        self.newickstring=''
        pass  