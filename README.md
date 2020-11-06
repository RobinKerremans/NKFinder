# NKFinder
Matlab software to derive optical constants from Transmission data

HOW TO:

1. Get R,T with the spectrophotometer. Make sure to be consistent with the 
position of the incoherent (Glass/Quartz) layer from the incoming light. Usually placed last.

2. Export the csv T/R data files to the NK_Finder map. Make sure they are in %

3. In the map, open NK_Finder.m

4. Only use the manual input region of the editor

5. Carefully read the commentary and initiate the variables properly

6. Select a rough range (ts ~50, thickness 30-600) 

7. Run the code. A graph will pop up to pinpoint the Cauchy regime

8. Click no and adjust the range of the cauchy regime if necessary

9. Run again. Click Yes this time.

10. A second graph will pop up showing the quality of the fit. Again there is the choice to update range.

11. Data of the rough fit will be displayed. Use this to run a finer fit (ts 5-10, narrower thickness range)

12. Continue and a csv file will be output with your specified name, containing n,k thickness and n_cauchy. 


