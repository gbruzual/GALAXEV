# GALAXEV: Galaxy Spectral Evolution Library
## (C) 2003-2023 - G. Bruzual and S. Charlot - All Rights Reserved
Currently, each **GALAXEV** SSP model is provided as a fits table file. Models are available for the values of the metallicity **Z** indicated in the table. The **2019** models (Plat et al. 2019; Sánchez et al. 2022) are based on the **PARSEC** evolutionary tracks (Chen et al. 2015; Marigo et al. 2013). In the **CB19** models (a.k.a. **C&B** models) the evolution of Post AGB stars follows Miller Bertolami (2018) whereas in the **BC19** models the Post AGB stars are treated as in Bruzual & Charlot (2003, **BC03**) following Vassiliadis & Wood (1994). The **BC03** models follow the **Padova 1994** evolutionary tracks and the **CB03** models the **Padova 2000** tracks (see BC03 for details). The **CB07** models are identical to the **BC03** models except for an enhanced contribution from TP-AGB stars (Bruzual 2007a,b).
<table>
<tr>
  <th>Model Release
  <th>2019
  <th>BC19/CB19
  <td>
  <th>Z:
  <td>0
  <td>
  <td>0.0001
  <td>
  <td>0.0002
  <td>
  <td>0.0005
  <th>
  <td>0.001
  <td>
  <td>0.002
  <td>
  <td>0.004
  <td>
  <td>0.006      
  <td>
  <td>0.008
  <td>
  <td>0.010
  <td>
  <td>0.014
  <td>
  <td>0.017
  <th>
  <td>0.02
  <td>
  <td>0.03
  <td>
  <td>0.04
  <td>
  <td>0.06
</tr>
<tr>
  <th>Model Release
  <th>2003
  <th>BC03/CB03/CB07
  <td>
  <th>Z:
  <td>
  <td>
  <td>0.0001
  <td>
  <td>
  <td>
  <td>0.0004
  <th>
  <td>0.001
  <td>
  <td>
  <td>
  <td>0.004
  <td>
  <td>     
  <td>
  <td>0.008
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <th>
  <td>0.02
  <td>
  <td>0.03
  <td>
  <td>0.05
  <td>
  <td>     
</tr>
</table>
Models have been computed for the combination of parameters listed in the following table. The <b>Stellar Library</b> entry refers to the stellar atlas used in the visible range, except for the BaSeL models which extend from the UV to the IR. To explore the <b>GALAXEV</b> models, several tasks are provided in the <b>pyGALAXEV</b> cell bellow. Default or recommended values are indicated in <span style="color:blue">blue</span>.
<table>
<tr>
  <th>Model
  <td>
  <th>IMF
  <td>
  <th>Mup
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <th>Model
  <td>
  <th>IMF
  <td>
  <th>Mup
  <td>
  <td>
  <td>
  <td>
  <td>
  <th>Stellar Libraries (all models)
<tr>
  <th><span style="color:blue">CB19</span>
  <td>
  <td>Chabrier, Kroupa, Salpeter
  <td>
  <td><span style="color:blue">100</span>, 10, 300, 600 Mo
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <th>BC19
  <td>
  <td>Chabrier
  <td>
  <td><span style="color:blue">100</span> Mo
  <td>
  <td>
  <td>
  <td>
  <td>
  <td><span style="color:blue">Miles</span>, Miles+, IndoUS, Stelib, BaSel
<tr>
  <th>
  <td>
  <td>Vazdekis et al. (1996)
  <td>
  <td><span style="color:blue">100</span> Mo
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <td>
  <th>BC03/CB03/CB07
  <td>
  <td>Chabrier, Kroupa, Salpeter
  <td>
  <td><span style="color:blue">100</span> Mo
  <td>
  <td>
<tr>
</table>
