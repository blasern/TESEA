context("triangle_score")

test_that("triangle_score runs", {
  if (requireNamespace("ESEA", quietly = TRUE)){
    # get data from ESEA
    ESEA::initializeESEA()
    edgesbackground <- ESEA::GetEdgesBackgrandData()
    pathwayEdge.db <- ESEA::GetPathwayEdgeData()
    dataset <- ESEA::GetExampleData("dataset")
    class.labels <- ESEA::GetExampleData("class.labels")
    controlcharacter <- ESEA::GetExampleData("controlcharactor") 
    # calculate edge score (triangle version)
    EdgeTriScore <- triangle_creation_score(dataset,
                                            class.labels, 
                                            controlcharacter, 
                                            edgesbackground[1:200, ])
    expect_true(EdgeTriScore['AANAT|ASMT'] == 0)
    expect_true(EdgeTriScore['AANAT|SLC25A16'] == 0)
    expect_error(triangle_creation_score(dataset,
                                         class.labels, 
                                         controlcharacter, 
                                         edgesbackground[1:10, ]), 
                 'Not enough overlap between genes in dataset and background network')
  }
})
