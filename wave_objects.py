from manimlib import * #type: ignore
import numpy as np


#functions
def get_light_color(wavelength, min= 1.9, max= 3.8):
    """
    Maps a wavelength (arbitrary units, e.g. ~200 nm per unit) to a Manim RGB color.
    Below min_wl → violet
    Above max_wl → red
    Between → smooth gradient.
    Returns a Color created from RGB values in [0, 1].
    """

    # Define the wavelength limits (in your arbitrary units)
    min_wl = 1.9
    max_wl = 3.5
    wl = np.clip(wavelength, min_wl, max_wl)#type:ignore

    # Convert to nanometers (roughly)
    wl_nm = np.interp(wl, [min_wl, max_wl], [380, 700])#type:ignore

    # --- Approximate visible spectrum RGB curve ---
    if 380 <= wl_nm < 440:
        R = -(wl_nm - 440) / (440 - 380)
        G = 0.0
        B = 1.0
    elif 440 <= wl_nm < 490:
        R = 0.0
        G = (wl_nm - 440) / (490 - 440)
        B = 1.0
    elif 490 <= wl_nm < 510:
        R = 0.0
        G = 1.0
        B = -(wl_nm - 510) / (510 - 490)
    elif 510 <= wl_nm < 580:
        R = (wl_nm - 510) / (580 - 510)
        G = 1.0
        B = 0.0
    elif 580 <= wl_nm < 645:
        R = 1.0
        G = -(wl_nm - 645) / (645 - 580)
        B = 0.0
    elif 645 <= wl_nm <= 700:
        R = 1.0
        G = 0.0
        B = 0.0
    else:
        R = G = B = 0.0

    # --- Apply simple intensity falloff near edges of visible spectrum ---
    if wl_nm > 700:
        factor = 0.3 + 0.7 * (780 - wl_nm) / (780 - 700)
    elif wl_nm < 380:
        factor = 0.3 + 0.7 * (wl_nm - 360) / (380 - 360)
    else:
        factor = 1.0

    R = np.clip(R * factor, 0, 1)#type:ignore
    G = np.clip(G * factor, 0, 1)#type:ignore
    B = np.clip(B * factor, 0, 1)#type:ignore

    # Return as Manim Color
    return Color(rgb=(R, G, B))

#test Scenes for testing new feaatures

"""class Test(Scene):
     def construct(self) -> None:
        axes = Axes(x_range=(-5,5,1),
                  y_range= (-5,5,1),
                  axis_config={
                      "include_numbers" : False
                  })
        def t_updater(mob, dt):
             mob.time += dt
        photon = WavePacket2D(axes= axes, wave_len= 800).add_updater(t_updater)
        self.add(photon)"""


class OscillatingWave2D(VMobject):
    def __init__(
            self,
            axes,
            amplitude = 1.0,
            wave_len = 2.0,
            speed = 1.0,
            phase = 0.0,
            wave_color= WHITE,
            arrow_density = 5,
            include_wave= True,
            include_arrows = True,
            arrow_opacity = .6,
            wave_opacity = 1.0,
            always_recalculate = False,
            wave_stroke_width=5,
            **kwargs
        ):
            """
            PARAMETERS
            =====================
            amplitude, wavelength, speed, phase, for basic wave build
            in the shape f(x,t) = a * sin(TAU/lambda * x - angular_freq * t - phase)
            (angular frequency calculated as: speed/wave_lenght * TAU)

            arrow_density: how many arrows in one unit,
            x_upper_bound : ValueTracker to manipulate the x_range upperbound to animate wave propogation
            x_lower_bound : animateable ValueTracker for the lowerbound of the x range

            always_recalculate: always recaltulate starting positions of vectors, useful if axes is continuisly moving in some way

            usage: manipulate time parameters with "wave".time and or x lower and upper bounds with "wave".x_lower_bound/upper_bound asd ValueTrackers
            """
            
            
            self.axes = axes
            self.amplitude = amplitude
            self.wave_len = wave_len
            self.speed = speed
            self.phase = phase
            self.wave_color = wave_color
            
            super().__init__(**kwargs)
            
            self.time = 0
            self.x_upper_bound = ValueTracker(self.axes.x_axis.x_max)
            self.x_lower_bound = ValueTracker(self.axes.x_axis.x_min)

            self.angular_freq = TAU * (self.speed / self.wave_len)
            
            self.wave_func = lambda x: self.amplitude * np.sin((TAU/self.wave_len) * x - self.angular_freq * self.time - self.phase)

            if include_wave:
                graph = always_redraw(
                    lambda:
                        self.axes.get_graph(self.wave_func, 
                                            stroke_color = self.wave_color, stroke_opacity = wave_opacity, stroke_width = wave_stroke_width,
                                            x_range= [self.x_lower_bound.get_value(), self.x_upper_bound.get_value()])
                )
            
                self.add(graph)#type: ignore
            
            max_arrow_num =  int((self.axes.x_axis.x_max - self.axes.x_axis.x_min) * arrow_density + 1)

            if include_arrows and always_recalculate:
                vects = always_redraw(lambda: VGroup(
                    Arrow(
                        start= self.axes.c2p(self.axes.x_axis.x_min, 0),
                        end= self.axes.c2p(self.axes.x_axis.x_min, 0) + self.axes.y_axis.get_unit_vector(),
                        buff = 0,
                        stroke_color = self.wave_color,
                        max_tip_length_to_length_ratio= .15
                    ).shift(self.axes.x_axis.get_unit_vector()* (1/arrow_density) * i).set_color(self.wave_color).set_opacity(arrow_opacity)
                    for i in range(max_arrow_num)
                ))

                def arrow_updater(vects):
                    for arrow in vects:
                        start = arrow.get_start()
                        if self.axes.p2c(start)[0] < self.x_upper_bound.get_value() and self.axes.p2c(start)[0] > self.x_lower_bound.get_value():
                            arrow.put_start_and_end_on(
                                start= start,
                                end= start + self.axes.y_axis.get_unit_vector() * self.amplitude * np.sin((TAU/self.wave_len) * self.axes.p2c(start)[0] - self.angular_freq * self.time - self.phase)
                            )
                        else:
                            arrow.put_start_and_end_on(
                                start= start,
                                end= start
                            )
                
                vects.add_updater(arrow_updater)

                self.add(vects)#type: ignore
            
            elif include_arrows and not always_recalculate:
                vects = VGroup(
                    Arrow(
                        start= self.axes.c2p(self.axes.x_axis.x_min, 0),
                        end= self.axes.c2p(self.axes.x_axis.x_min, 0) + self.axes.y_axis.get_unit_vector(),
                        buff = 0,
                        stroke_color = self.wave_color,
                        max_tip_length_to_length_ratio= .15
                    ).shift(self.axes.x_axis.get_unit_vector()* (1/arrow_density) * i).set_color(self.wave_color).set_opacity(arrow_opacity)
                    for i in range(max_arrow_num)
                    )

                def arrow_updater(vects):
                    for arrow in vects:
                        start = arrow.get_start()
                        if self.axes.p2c(start)[0] < self.x_upper_bound.get_value() and self.axes.p2c(start)[0] > self.x_lower_bound.get_value():
                            arrow.put_start_and_end_on(
                                start= start,
                                end= start + self.axes.y_axis.get_unit_vector() * self.amplitude * np.sin((TAU/self.wave_len) * self.axes.p2c(start)[0] - self.angular_freq * self.time - self.phase)
                            )
                        else:
                            arrow.put_start_and_end_on(
                                start= start,
                                end= start
                            )
                
                vects.add_updater(arrow_updater)
                self.add(vects)

        

class WavePacket2D(OscillatingWave2D):
    def __init__(
            self, 
            axes, 
            amplitude=1.0, 
            wave_len=400, 
            speed=1.0, 
            phase=0.0, 
            wave_color=WHITE, 
            arrow_density=10, 
            include_wave=True, 
            include_arrows=False, 
            arrow_opacity=0.6, 
            wave_opacity=1.0,
            always_recalculate = False, 
            **kwargs
            ):
                super().__init__(axes, amplitude, wave_len/400, speed,
                                phase, get_light_color(wavelength= wave_len/200), arrow_density, include_wave,
                                include_arrows, arrow_opacity, wave_opacity, always_recalculate, 8, **kwargs)

                self.time = 0
                
                
                
                hull = Ellipse(
                        width= self.amplitude * 4.7, 
                        height= self.amplitude * 3.7, 
                        stroke_color = self.wave_color,
                        stroke_width = 10 * self.amplitude, 
                        fill_color = self.wave_color, 
                        fill_opacity = .3
                        )
                
                def boundry_upd(mob):
                    mob.x_lower_bound.set_value(self.axes.p2c(hull.get_edge_center(self.axes.x_axis.get_unit_vector() * (-1)))[0] + (0.5) )
                    mob.x_upper_bound.set_value(self.axes.p2c(hull.get_edge_center(self.axes.x_axis.get_unit_vector()))[0] - (0.5)) 
                
                self.add_updater(boundry_upd)
                
                self.add(hull)
          
class WaveSum2D(VMobject):
    def __init__(
            self,
            axes,
            waves_list,
            wave_color=WHITE,
            arrow_density = 5,
            include_wave= True,
            include_arrows = True,
            arrow_opacity = .6,
            wave_opacity = 1.0,
            always_recalculate = False,
            **kwargs    
    ):
        super().__init__(**kwargs)
        self.axes = axes
        self.waves_list = waves_list
        self.wave_color = wave_color
        
        self.x_upper_bound = ValueTracker(self.axes.x_axis.x_max)
        self.x_lower_bound = ValueTracker(self.axes.x_axis.x_min)

        def sum_func(x):
            y = 0 
            for wave in waves_list:
                y += wave.wave_func(x)
            return y    
        
        if include_wave:
            graph = always_redraw(lambda: self.axes.get_graph(sum_func,
                                                            stroke_color = self.wave_color, stroke_opacity = wave_opacity,
                                                            x_range= [self.x_lower_bound.get_value(), self.x_upper_bound.get_value()]))
            self.add(graph)#type: ignore

        max_arrow_num =  (self.axes.x_axis.x_max - self.axes.x_axis.x_min) * arrow_density + 1

        if include_arrows and always_recalculate:
            vects = always_redraw(lambda: VGroup(
                Arrow(
                    start= self.axes.c2p(self.axes.x_axis.x_min, 0),
                    end= self.axes.c2p(self.axes.x_axis.x_min, 0) + self.axes.y_axis.get_unit_vector(),
                    buff = 0,
                    stroke_color = self.wave_color,
                    max_tip_length_to_length_ratio= .15
                    ).shift(self.axes.x_axis.get_unit_vector()* (1/arrow_density) * i).set_color(self.wave_color).set_opacity(arrow_opacity)
                    for i in range(max_arrow_num)
                ))

            def arrow_updater(vects):
                for arrow in vects:
                    start = arrow.get_start()
                    if self.axes.p2c(start)[0] < self.x_upper_bound.get_value() and self.axes.p2c(start)[0] > self.x_lower_bound.get_value():
                        arrow.put_start_and_end_on(
                            start= start,
                            end= start + self.axes.y_axis.get_unit_vector() * sum_func(self.axes.p2c(start)[0])
                        )
                    else:
                        arrow.put_start_and_end_on(
                            start= start,
                            end= start
                            )
                
                vects.add_updater(arrow_updater)

                self.add(vects)#type: ignore
            
        elif include_arrows and not always_recalculate:
            vects = VGroup(
                    Arrow(
                        start= self.axes.c2p(self.axes.x_axis.x_min, 0),
                        end= self.axes.c2p(self.axes.x_axis.x_min, 0) + self.axes.y_axis.get_unit_vector(),
                        buff = 0,
                        stroke_color = self.wave_color,
                        max_tip_length_to_length_ratio= .15
                    ).shift(self.axes.x_axis.get_unit_vector()* (1/arrow_density) * i).set_color(self.wave_color).set_opacity(arrow_opacity)
                    for i in range(max_arrow_num)
                    )

            def arrow_updater(vects):
                for arrow in vects:
                    start = arrow.get_start()
                    if self.axes.p2c(start)[0] < self.x_upper_bound.get_value() and self.axes.p2c(start)[0] > self.x_lower_bound.get_value():
                        arrow.put_start_and_end_on(
                            start= start,
                            end= start + self.axes.y_axis.get_unit_vector() * sum_func(self.axes.p2c(start)[0])
                        )
                    else:
                        arrow.put_start_and_end_on(
                            start= start,
                            end= start
                        )
                
            vects.add_updater(arrow_updater)
            self.add(vects)
