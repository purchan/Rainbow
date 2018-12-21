# Rainbow
## Target
- Multiple Drop
- Camera Spin
- Multiple "Bag of Light"
  - Light : Color and Refraction Index Correction
- Physics
  - Refraction
  - Reflection
  
## Reference
### Light : Color and Refraction Index Correction
- [Refractive Equations](https://en.wikipedia.org/wiki/Sellmeier_equation)
- [RefractiveIndex.INFO](https://refractiveindex.info)
  - [Air](https://refractiveindex.info/?shelf=other&book=air&page=Ciddor)
    - P. E. Ciddor. Refractive index of air: new equations for the visible and near infrared, [Appl. Optics 35, 1566-1573 (1996)](https://doi.org/10.1364/AO.35.001566)
    - [Calculation script (Python)](https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Ciddor%201996%20-%20air.py) - can be used for calculating refractive index of air at a given humidity, temperatire, pressure, and CO2 concentration
  - [Water](https://refractiveindex.info/?shelf=main&book=H2O&page=Daimon-19.0C)
    - M. Daimon and A. Masumura. Measurement of the refractive index of distilled water from the near-infrared region to the ultraviolet region, [Appl. Opt. 46, 3811-3820 (2007)](https://doi.org/10.1364/AO.46.003811)
- [Light Wavelength to RGB](https://academo.org/demos/wavelength-to-colour-relationship/)
