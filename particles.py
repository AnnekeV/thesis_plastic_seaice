
class Sea_Ice_Particle(JITParticle):         # Define a new particle class to sample sea ice
    sit        = Variable('sit', initial=fieldset.sit)    # Sea ice thickness
    sic        = Variable('sic', initial=fieldset.sic)    # Sea ice concentration
    in_ice     = Variable('in_ice', initial = 1.)
    prev_sit   = Variable('prev_sit', initial=fieldset.sit)
    prev_state = Variable('prev_state', initial=1.)
