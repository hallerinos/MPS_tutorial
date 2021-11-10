# Welcome!

Before you start digging in the exercises, let's fix some necessities right away.
First of all, we must choose a unique ordering of our tensors that we want to contract in a network.
The choice is fully up to you, but in the solutions you will find the following one:

```
 ----                               ---- 
|    |--(-1)                 (-1)--|    |
|    |                             |    |
|    |                             |    |
|    |                             |    |
|    |                             |    |
|    |                             |    |
| LW |--(-2)                 (-2)--| RW |
|    |                             |    |
|    |                             |    |
|    |                             |    |
|    |                             |    |
|    |                             |    |
|    |--(-3)                 (-3)--|    |
 ----                               ---- 

                                    (-4)      
                                     |        
       -----------                  ---       
(-1)--|    AC     |--(-4)    (-1)--| W |--(-3)
       -----------                  ---       
        |       |                    |        
       (-2)    (-3)                 (-2)      
```