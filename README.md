# FlyLimbOptoCode

Code to perform localized optogenetic stimulation using a simple method, as described in our [paper](https://doi.org/10.7554/eLife.54183). If you find this helpful, please consider citing: DeAngelis, B.D.\*, Zavatone-Veth, J.A.\*, Gonzalez-Suarez, A.D., and Clark, D.A. Spatiotemporally precise optogenetic activation of sensory neurons in freely walking *Drosophila*. *eLife* 9: e54183 (2020). DOI: [10.7554/eLife.54183](https://doi.org/10.7554/eLife.54183). (\*equal contributions).

## Dependencies

The code to generate and project patterned stimuli relies on several packages and libraries, all of which are freely available.

- [Psychtoolbox](http://psychtoolbox.org/) is used for stimulus generation and display control. 
- Bindings for the [LightCrafter 4500](http://www.ti.com/tool/DLPLCR4500EVM) projector are provided by [matlab-lcr](https://github.com/Stage-VSS/matlab-lcr), an [MIT-licensed](https://opensource.org/licenses/MIT) library. For convenience, we include a copy under [Projection/externals/matlab-lcr](https://github.com/ClarkLabCode/FlyLimbOptoCode/tree/master/Projection/externals/matlab-lcr).
- Bindings for other LightCrafter projectors are provided by the [TI-DLP-LightCrafter](https://github.com/fglichttechnik/TI-DLP-LightCrafter) library, which is [free software](https://github.com/fglichttechnik/TI-DLP-LightCrafter/blob/master/README.txt). For convenience, we include a copy of the relevant code under [Projection/externals](https://github.com/ClarkLabCode/FlyLimbOptoCode/tree/master/Projection/externals).


Brian D DeAngelis1 et al (2020). Spatiotemporally precise optogenetic activation of sensory neurons in freely walking Drosophila. DOI: 10.7554/eLife.54183.
