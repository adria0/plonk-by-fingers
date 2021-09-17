# plonk-by-fingers

This is a quick & dirty implementation of the excellent Joshua Fitzgerald [Plonk by hand](https://research.metastate.dev/plonk-by-hand-part-1) ([part2](https://research.metastate.dev/plonk-by-hand-part-2-the-proof)) ([part3](https://research.metastate.dev/plonk-by-hand-part-3-verification)) tutorial

- do not expect this code to be anything close to production, is intended just to understand the protocol
- there is a mistake in the hand computations in part3 that says that <img src="https://render.githubusercontent.com/render/math?math=l_{16P,P}=x%2B15"> where the correct value seems to be <img src="https://render.githubusercontent.com/render/math?math=l_{16P,P}=x%2B100">, this also affects the pairing, that is <img src="https://render.githubusercontent.com/render/math?math=93%2b76u"> instead <img src="https://render.githubusercontent.com/render/math?math=97%2b89u"> 
