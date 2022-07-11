from manim import *

class displayEquations(Scene):
    def construct(self):
        # Create Tex objects
        eq1 = ['$Y_i =$', '$\\beta_1 X_i + \\beta_2 D_i +$', '$ \\epsilon_i$']
        eq2 = ['$Y_i  =$', '$\\beta_1 X_i + \\beta_2 D_i +$', '$ U_i$','$+$','$\\epsilon_i$']
        eq3 = ['$Y_{it} =$', '$\\beta_1 X_{it} +$', '$\\beta_2 D_{it}+$', '$U_{i}$', '+','$\\epsilon_{it}$']
        eq4 = ['$Y_{i,t-1} =$', '$\\beta_1 X_{i,t-1} +$','$\\beta_2 D_{i,t-1}+$', '$U_{i}$','$+\\epsilon_{i,t-1}$']
        eq5 = ['$(Y_{it} - Y_{i,t-1}) =$', '$\\beta_1 (X_{it}-X_{i,t-1}) +$','$\\beta_2 (D_{it}-D_{i,t-1})+$', '$(U_{i}-U_i)$','$+ (\\epsilon_{it}-\\epsilon_{i,t-1})$']
        eq6 = ['$(Y_{it} - Y_{i,t-1}) =$', '$\\beta_1 (X_{it}-X_{i,t-1}) +$','$\\beta_2 (D_{it}-D_{i,t-1})$','+','$ (\\epsilon_{it}-\\epsilon_{i,t-1})$']
        eq7 = ['$\\Delta Y_{it} =$', '$\\beta_1 \\Delta X_{it} +$', '$\\beta_2 \\Delta D_{it}+$','$\\Delta \\epsilon_{it}$']

        first_line = Text("Basic regression")
        second_line = Text("Add individual heterogeneity...")
        third_line = Text("Move to panel data setting ...")
        fourth_line = Text("Define for two time periods...")
        fifth_line = Text("Take the difference...")
        
        t1 = VGroup(*[Tex(eq) for eq in eq1]).arrange(RIGHT, buff=0.2)
        t2 = VGroup(*[Tex(eq) for eq in eq2]).arrange(RIGHT, buff=0.2)
        t3 = VGroup(*[Tex(eq) for eq in eq3]).arrange(RIGHT, buff=0.2)
        t4 = VGroup(*[Tex(eq) for eq in eq4]).arrange(RIGHT, buff=0.2)
        t5 = VGroup(*[Tex(eq) for eq in eq5]).arrange(RIGHT, buff=0.2)
        t6 = VGroup(*[Tex(eq) for eq in eq6]).arrange(RIGHT, buff=0.2)
        t7 = VGroup(*[Tex(eq) for eq in eq7]).arrange(RIGHT, buff=0.2)
        
        # Position second line underneath first line
        t3.next_to(t2, DOWN)
        t4.next_to(t3, DOWN)
        t5.next_to(t2,DOWN)
        t6.next_to(t2,DOWN)
        t7.next_to(t3,DOWN)
        
        first_line.next_to(t1,UP)
        second_line.next_to(t1,UP)
        third_line.next_to(t3,DOWN)
        fourth_line.next_to(t1,UP)
        fifth_line.next_to(t1,UP)
        
        t2[2].set_color(RED) 
        t3[3].set_color(RED) 
        t4[3].set_color(RED) 
        t5.scale(0.8)
        t6.scale(0.8)
        
        # Displaying text and equation
        # Basic regression
        self.play(Write(t1),Write(first_line))
        self.wait(1)
        
        # Add individual heterogeneity
        self.play(Transform(first_line,second_line))
        self.wait(1)
        self.play(TransformFromCopy(t1[0:2],t2[0:2]),TransformFromCopy(t1[2], t2[2:5]),FadeOut(t1))
        self.wait(1)
        
        # Move to panel data setting
        self.play(FadeIn(third_line))
        self.play(TransformFromCopy(t2,t3))
        self.wait(2)
        self.play(FadeOut(third_line))
        self.play(FadeOut(t2),FadeOut(second_line))
        self.play(TransformFromCopy(t3,t4))
        self.wait(2)
        self.play(Transform(t3[0],t5[0]),FadeOut(t4[0]))
        self.play(ReplacementTransform(t4[2],t5[1]),FadeOut(t3[1]))
        self.play(ReplacementTransform(t4[3],t5[2]),FadeOut(t3[2]),ReplacementTransform(t4[4],t5[3]),FadeOut(t3[3]),ReplacementTransform(t4[5],t5[4]),FadeOut(t3[4],t3[5]))
        self.play(FadeOut(t5[3]))
        self.wait(1)
        
        self.play(TransformFromCopy(t5[0],t7[0]))
        self.play(TransformFromCopy(t5[1],t7[1]))
        self.play(TransformFromCopy(t5[2],t7[2]))
        self.play(TransformFromCopy(t5[4],t7[3]))
        
        
        

        #self.play(ReplacementTransform(t2,t3))
        self.wait(3)
