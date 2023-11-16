# Models
This folder contains MATLAB models for image classification. The supported models include 'TreeBagger', 'ClassificationKNN', and 'ClassificationECOC'. Additionally, you can use a custom MATLAB struct adhering to the following format:

```model.probAllImage: % Probabilities
model.classMap: % Predicted labels
model.listOfClass: % Probabilities order```

## Supported Models

1. TreeBagger [using TreeBagger](https://www.mathworks.com/help/stats/treebagger.html)
2. ClassificationKNN  [using fitcknn](https://www.mathworks.com/help/stats/fitcknn.html)
3. ClassificationECOC [using fitcsvm](https://www.mathworks.com/help/stats/fitcecoc.html)

## Custom Model Struct

If you choose to use a custom model struct, ensure that it includes the following fields:

* model.probAllImage: This field should contain the probabilities associated with each class for each image.
* model.classMap: Predicted labels for each image.
* model.listOfClass: The order of probabilities corresponding to each class.
