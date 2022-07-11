from manim import *

class ExampleEquation(Scene): 
    def construct(self):
        self.method_a()
        # self.method_b()                                                                                                                                                                                                                      

    def method_a(self):
        eq1 = ["$Y_i$ =", "$X_i'\\beta + \\tau D_i$", "$+ \\epsilon_i$"]
        eq2 = ["$Y_i$ =", "$X_i'\\beta + \\tau D_i$", "$+ U_i + \\epsilon_i$"]
        eq3 = ["$Y_i$ =", "$X_i'\\beta + \\tau D_i$", "$+ \\underbrace{U_i}_{\\text{Heterogeneity}} + \\epsilon_i$"]        

        t1 = VGroup(*[Tex(eq) for eq in eq1]).arrange(RIGHT, buff=0.2)

        t2 = VGroup(*[Tex(eq) for eq in eq2]).arrange(RIGHT, buff=0.2)

        t3 = VGroup(*[Tex(eq) for eq in eq3]).arrange(RIGHT, buff=0.2)        

        t2.next_to(t1, DOWN)
        t3.next_to(t2, DOWN)        

#        t1[2].set_color(RED) #b
#        t2[2].set_color(RED) #c+d
#        t2[2].set_color(RED)
        
        self.play(Write(t1), run_time=2)
        self.play(LaggedStart(
            TransformFromCopy(t1[0], t2[0]),
            TransformFromCopy(t1[1], t2[1]),
            TransformFromCopy(t1[2], t2[2]),
            TransformFromCopy(t1[2], t3[2]),            
            lag_ratio=0.9,
            run_time=5,
        ))
        self.wait()
