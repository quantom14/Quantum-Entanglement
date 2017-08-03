Examples 
===

Here are a few lines of code that demonstrate how NegMapEns works. We create to cells, rho and sigma and we fill them with product states and maximally entangled states. At each step we compute how much Negativity we need if we wish to transform the input states into the output states.

```
>> rho = cell(1,4);
>> sigma = cell(1,4);
>> k = 1;
>> for i = 0:1
    	for j = 0:1
        	rho{k} = prod_state(i,j)*prod_state(i,j)';
        	sigma{k} = Bell(k-1)*Bell(k-1)';
        
        	NegMapEns(rho(1:k),2,sigma(1:k),2)
        
       		k = k+1; 
    	end
   end

ans =

    0.1250


ans =

    0.2500


ans =

    0.3750


ans =

    0.5000
````
If we relax the mapping to a non-zero value of epsilon 
````
>> NegMapEns(rho,2,sigma,2,0.5)

ans =

    0.2500

````