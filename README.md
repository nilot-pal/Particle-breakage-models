This repository contains supporting codes and slides for particle breakage stochastic models. 
# The problem
Aero engine manufacturing companies do extensive CFD analysis of air flow inside an engine to design more effcient air-breathing engines. However, all this flow analysis assumes air to be clean and devoid of solid particles. This is done to avoid all the complexities of fluid-particle, particle-particle coupling in CFD. 
<p align="center">
<img src="https://user-images.githubusercontent.com/72824334/209989558-2bbd9fee-d2ba-4499-9e6b-c9a7e59246a1.png" width="400" height="400">
</p>  
<p align="center">
    <em>Ice crystals mixed with incoming air enter into the engine</em>
</p>

However in cold regions such as Finland, before an aircraft is taking off, particles are sprinkled on the runway to melt the ice. When the engine is started and the aircraft is moving on the runway, these particles are sucked into the engine and can impact several regions inside. Particles can break into smaller fragments on impact
with the compressor blades (say) and cause compressor surge or stall. Similarly, particles can deposit on blades and if this happens over a long time, engine gets damaged (Although cleaning particle deposits using ultrasonic bath is an option, but nevertheless it can't remove all particles). In the hotter middle eastern regions, there is no snow but the atmospheric air contains a large concentration of sand and dust. While cruising, an aircraft easily ingests these particles which subsequently damage the engines. Even temperate regions like India pose a threat to these engines due to high concentration of pollutants in the atmosphere (New Delhi tops the list).

