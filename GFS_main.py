from manimlib import * #type: ignore
import numpy as np
import sys, os

sys.path.append(os.path.dirname("C:\\Users\\Priyam\\ManimGL\\manim\\Own\\GFS"))

from GFS.photoelectric_effect import *
from GFS.wave_objects import *

class MacroscopicExperiment(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        electroscope  = Electroscope(scale= 1.5).shift(DOWN + 3* RIGHT)#type: ignore
        light_source = LightCone(wave_length=400,auto_create=False)
        light_source.source.shift(1.15*UP + 3 * RIGHT)#type: ignore
        light_intensity = ValueTracker(0)
        light = always_redraw(lambda : Polygon(light_source.source[0].get_center(), [3,1.65,0], [3,1.1,0], stroke_opacity = 0, fill_color = light_source.get_light_color(light_source.wave_len), fill_opacity = light_intensity.get_value() / 2))#type: ignore
        
        self.add(light, electroscope ,light_source.source)
         
        #self.wait(3)
        self.play(electroscope.charge.animate(rate_func = overshoot, run_time = 3).set_value(-0.7))
        #self.wait()
        #light_intensity.set_value(0.7)
        #self.play(electroscope.charge.animate(rate_func = overshoot, run_time = 3).set_value(0))
        #light_intensity.set_value(0)
        #self.wait(2)

class SimpleExperimentStills(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        experiment = PhotoElectricExperimentBasic(incoming_light_wave_len= 420,
                                                  incoming_light_intensity= 0.0,
                                                  electron_lifetime= 14,
                                                  electron_x_boundary=(-7,4),
                                                  performant_eField= True,
                                                  always_redraw_e_field= False,
                                                  emitter_metal= "Cs",
                                                  simulation_speed= 1,)
        experiment.collector_top.shift(6*LEFT)#type: ignore
        experiment.collector_bottom.shift(6*LEFT)#type: ignore
        
        self.add(experiment)
        

class SimpleExperimentSpeedVects(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        experiment = PhotoElectricExperimentBasic(incoming_light_wave_len= 420,
                                                  incoming_light_intensity= 0.5,
                                                  electron_lifetime= 14,
                                                  electron_x_boundary=(-7,4),
                                                  performant_eField= True,
                                                  always_redraw_e_field= False,
                                                  emitter_metal= "Cs",
                                                  simulation_speed= 1,)
        experiment.collector_top.shift(6*LEFT)#type: ignore
        experiment.collector_bottom.shift(6*LEFT)#type: ignore
        
        
        lightintensity = ValueTracker(0)

        def intensity_upd(emitter):
            emitter.emission_intensity = lightintensity.get_value()
        def light_op_upd(light):
            light.intensity = lightintensity.get_value()
        
        experiment.emitter.add_updater(intensity_upd)
        experiment.lightsource.add_updater(light_op_upd)
        self.add(experiment)
        
        #self.wait(1)
        lightintensity.set_value(.0)
        light_on = lightintensity.animate(rate_func = linear, run_time = 0.5).set_value(0.5)
        pause_8 = Pause(experiment, run_time= 3)
        
        self.play(LaggedStart(light_on, pause_8, lag_ratio= 1))
        experiment.emitter.clear_updaters()
        
        vect_anim_group = []
        
        for electron in experiment.emitter.electron_group:
            
            vect_buff = 0.065 * normalize(electron.velocity)#type: ignore
            v_vect = Arrow(start= electron.get_center() + vect_buff, 
                           end= electron.get_center() + vect_buff + (electron.velocity * 0.2), #type: ignore
                           buff= 0,
                           max_width_to_length_ratio= 0.05,
                           max_tip_length_to_length_ratio= 0.15,
                           ).set_color(GREY_A).set_opacity(0.8)
            vect_anim_group.append(GrowArrow(v_vect))
            
        
        #self.play(self.camera.frame.animate.shift(RIGHT), run_time= 0.5)
        #self.play(self.camera.frame.animate.scale(0.25), run_time = 1)
        #self.wait(2)
        
        self.play(*vect_anim_group, run_time = 1)
        #self.play(light_on, pause_8)
        #lightintensity.set_value(0)
        #self.wait(10)


class Test(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        

class Pause(Animation):
    def __init__(self, 
                 mobject, 
                 run_time= 1, 
                 time_span= None, 
                 lag_ratio= 0, 
                 rate_func= linear, 
                 name = "", 
                 remover = False, 
                 final_alpha_value = 1, 
                 suspend_mobject_updating = False
                 ):
        super().__init__(mobject, run_time, time_span, lag_ratio, rate_func, name, remover, final_alpha_value, suspend_mobject_updating)     
              

#WaveModelContradiction
class WaveModelVis(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        

        def t_upd(mob, dt):
            mob.time += dt

        wave_lmd = 2 # im 200 nm
        wave_speed = 4
        wave_amplitude = 2.4 # if amplitude of the wave is above 2, then the electron will get ejected
        wave_phase = 1*PI
        
        incoming_light_angle = 10 * DEG
        ejection = False
        ejection_time = 8.5/(wave_speed/3)
        sim_time_step = 0.15 * ( wave_speed/3)* (4/wave_lmd) # for every one unit wavespeed add 0.05 to the timestep

        if wave_amplitude > 2:
            ejection = True

        self.camera.frame.scale(1.3).shift(0*LEFT)#type: ignore
        #fitrst scene shows a wave with a low amplitude unsuccesfully causing an electron to wiggle, but not realease
        wave_axes = Axes(x_range=(-7,11,1),
                         y_range=(-5,5,1)).rotate(angle= incoming_light_angle, about_point= (5,0,0), axis= Z_AXIS)#type:ignore
        wave_low_a = OscillatingWave2D(axes=wave_axes,
                                       amplitude=wave_amplitude, 
                                       wave_color= get_light_color(wave_lmd), 
                                       wave_len= wave_lmd, 
                                       speed= wave_speed, 
                                       phase= wave_phase,
                                       include_arrows= True,
                                       ).add_updater(t_upd)
        
        metal = Rectangle(width= 3, height= 20, stroke_opacity = 0, fill_color = GREY_D, fill_opacity = 1).shift(RIGHT *6)#type: ignore

        equillibrium_line = DashedVMobject(wave_axes.get_graph(lambda x: 0, x_range=[-7,20],stroke_opacity = 0.5, stroke_color = GREY_A), num_dashes= 100)
        wave_low_a.x_upper_bound.set_value(-6)

        def e_pos_upd_low(e,dt):
            e.time += dt
            if wave_axes.p2c(e.get_center())[0] < wave_low_a.x_upper_bound.get_value() and not ejection:#type: ignore
                e.acceleration = wave_axes.y_axis.get_unit_vector()*wave_low_a.wave_func(wave_axes.p2c(e.get_center())[0])
                e.shift(e.acceleration  * sim_time_step)
            elif wave_axes.p2c(e.get_center())[0] < wave_low_a.x_upper_bound.get_value() and ejection:#type:ignore
                if e.time < ejection_time:
                    e.acceleration = wave_axes.y_axis.get_unit_vector()*wave_low_a.wave_func(wave_axes.p2c(e.get_center())[0])
                    e.shift(e.acceleration  * sim_time_step)
                else:
                    e.acceleration += sim_time_step * 0.25 * LEFT#type:ignore
                    e.shift(e.acceleration  * sim_time_step)
        
        def a_vect_upd(vect):
            if not ejection:
                if wave_axes.p2c(electron.get_center())[0] < wave_low_a.x_upper_bound.get_value():#type: ignore
                    vect.put_start_and_end_on(start = electron.get_center(), end= electron.get_center() + wave_axes.y_axis.get_unit_vector() * wave_low_a.wave_func(wave_axes.p2c(electron.get_center())[0]))
                else:
                    vect.put_start_and_end_on(ORIGIN, ORIGIN)
            elif ejection:
                if electron.time < ejection_time:
                    if wave_axes.p2c(electron.get_center())[0] < wave_low_a.x_upper_bound.get_value():#type: ignore
                        vect.put_start_and_end_on(start = electron.get_center(), end= electron.get_center() + wave_axes.y_axis.get_unit_vector() * wave_low_a.wave_func(wave_axes.p2c(electron.get_center())[0]))
                    else:
                        vect.put_start_and_end_on(ORIGIN, ORIGIN)
                else:
                    vect.put_start_and_end_on(start= electron.get_center(), end = electron.get_center() + electron.acceleration)

        electron = Electron().shift(5*RIGHT-wave_axes.y_axis.get_unit_vector()*wave_low_a.amplitude).add_updater(e_pos_upd_low)#type:ignore
        acceleration_vect = Arrow(start=ORIGIN, end= UP, buff = 0).set_color(GREEN).add_updater(a_vect_upd)


        self.add(metal, wave_low_a, electron, acceleration_vect, equillibrium_line)
        self.play(wave_low_a.x_upper_bound.animate.set_value(11),run_time = 17/wave_speed, rate_func = linear)
        self.wait(10)

        
        

class BasicSimIntro(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        wave_len = 400
        intensity = 0.0
        final_voltage = 1.16
        emitter_metal = "Cs"
        simulation_speed = 1


        experiment_basic = PhotoElectricExperimentBasic(incoming_light_intensity= intensity, 
                                                        incoming_light_wave_len= wave_len,
                                                        emitter_metal= emitter_metal, 
                                                        performant_eField=True, 
                                                        simulation_speed= simulation_speed,
                                                        always_redraw_e_field=False,
                                                        )
        
        electron_speed_label = DecimalNumber(number= experiment_basic.electron_speed_real, num_decimal_places=0).to_corner(UL)
        m_s_label = Text("m/s").to_corner(UL).shift(2.7 * RIGHT + 0.05 * UP)#type: ignore
        voltage_label = always_redraw(lambda: DecimalNumber(number= experiment_basic.voltage.get_value(), num_decimal_places= 2).to_corner(UR).shift(0.5*LEFT))#type: ignore
        V_label = Text("V").to_corner(UR)
        wavelen_label = DecimalNumber(number= experiment_basic.incoming_light_wave_len, num_decimal_places= 1).to_edge(UP)

        #self.add(electron_speed_label, m_s_label, voltage_label, V_label, wavelen_label)
        self.add( experiment_basic)
        #self.wait(3)
        #self.play(experiment_basic.voltage.animate(run_time = 6).set_value(final_voltage)) 
      
class BasicSimulation(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = [0.05,0.05,0.05,1]
        
        wave_len = 500
        intensity = ValueTracker(0.0)
        final_intensity = 0.65
        final_voltage = 0.537#3.02#250 #1.165# 400 # 0.537#500
        emitter_metal = "Cs"
        simulation_speed = 1

        

        experiment_basic = PhotoElectricExperimentBasic(
                                                        incoming_light_intensity= intensity.get_value(), #type:ignore
                                                        incoming_light_wave_len= wave_len,
                                                        emitter_metal= emitter_metal, 
                                                        performant_eField=True, 
                                                        simulation_speed= simulation_speed,
                                                        always_redraw_e_field=False,
                                                        )
        
        e_field_show  = VGroup(
                        Arrow(start= (1.9, 1.5 - 0.5*i, 0),#type:ignore
                              end = (-1.9, 1.5 - 0.5*i, 0),#type:ignore
                              buff = 0,
                              ).set_color(GREEN_D)
                            for i in range(7)
                        )
        
        def field_upd(e_field):
            for arr in e_field:
                arr.set_opacity((experiment_basic.voltage.get_value() / final_voltage) * 0.25)
        e_field_show.add_updater(field_upd)
        
        
        def intensity_upd(emitter):
            emitter.emission_intensity = intensity.get_value()
        def light_op_upd(light):
            light.intensity = intensity.get_value()
        
        experiment_basic.emitter.add_updater(intensity_upd)
        experiment_basic.lightsource.add_updater(light_op_upd)

        current_meter = Amperemeter(init_current= intensity.get_value(), scale= 0.35).move_to((-2, -3, 0)) #type: ignore
        
        electron_speed_label = DecimalNumber(number= experiment_basic.electron_speed_real, num_decimal_places=0).to_corner(UL)
        m_s_label = Text("m/s").to_corner(UL).shift(2.7 * RIGHT + 0.05 * UP)#type: ignore
        voltage_label = always_redraw(lambda: DecimalNumber(number= experiment_basic.voltage.get_value(), num_decimal_places= 2).to_corner(UR).shift(LEFT + 2* DOWN))#type: ignore
        V_label = Text("V").to_corner(UR)
        wavelen_label = DecimalNumber(number= experiment_basic.incoming_light_wave_len, num_decimal_places= 0).move_to((-5, 1, 0))#type:ignore

        #self.add(electron_speed_label, m_s_label, voltage_label, V_label, wavelen_label)
        self.add( experiment_basic, current_meter, voltage_label,wavelen_label, e_field_show)
        light_anim = intensity.animate(run_time = 0.5, rate_func = linear).set_value(final_intensity)
        meter_anim_pos = current_meter.current.animate(run_time= final_intensity * 8, rate_func = overshoot).set_value(final_intensity)
        #meter_anim_neg = current_meter.current.animate(run_tim = final_intensity * 10, rate_func = overshoot).set_value(0)
        
        self.play(light_anim)
        self.wait(2)
        self.play(meter_anim_pos)
        self.wait(3)
        self.play(LaggedStart(experiment_basic.voltage.animate(run_time = 4, rate_func = linear).set_value(final_voltage),
                  current_meter.current.animate(run_time = final_intensity * 10, rate_func = smooth).set_value(0),
                  lag_ratio=0.4))
        self.wait(10)
        
class AccurateSimulation(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = (0.05,0.05,0.05,1)     #type: ignore
        

        wave_len = 700
        intensity = ValueTracker(0.0)
        final_intensity = .6
        final_voltage = 0.0#3.02#250 #1.165# 400 # 0.537#500
        emitter_metal = "Cs"
        simulation_speed = 1

        experiment = PhotoElectricExperimentPhoton(incoming_light_wave_len=wave_len, 
                                                   incoming_light_intensity= 0, 
                                                   performant_eField= False,
                                                   electron_lifetime= 13,
                                                   simulation_speed= 1,
                                                   emitter_metal= "Cs"
                                                   )
        
        current_meter = Amperemeter(init_current= intensity.get_value(), scale= 0.35).move_to((-2, -3, 0)) #type: ignore

        

        def intensity_upd(experiment_):
            experiment_.light_source.emission_intensity = intensity.get_value()
        experiment.add_updater(intensity_upd)
        
        e_field_show  = VGroup(
                        Arrow(start= (1.9, 1.5 - 0.5*i, 0),#type:ignore
                              end = (-1.9, 1.5 - 0.5*i, 0),#type:ignore
                              buff = 0,
                              ).set_color(GREEN_D)
                            for i in range(7)
                        )
        
        def field_upd(e_field):
            for arr in e_field:
                arr.set_opacity((experiment.voltage.get_value() / final_voltage) * 0.25)
        e_field_show.add_updater(field_upd)

        voltage_label = always_redraw(lambda: DecimalNumber(number= experiment.voltage.get_value(), num_decimal_places= 2).to_corner(UR).shift(LEFT + 2* DOWN))#type: ignore
        wavelen_label = DecimalNumber(number= experiment.incoming_light_wave_len, num_decimal_places= 0).move_to((-5, 1, 0))#type:ignore

        light_anim = intensity.animate(run_time = 0.5, rate_func = linear).set_value(final_intensity)
        meter_anim_pos = current_meter.current.animate(run_time= final_intensity * 6, rate_func = overshoot).set_value(final_intensity * 0.5)


        self.add(experiment, e_field_show, voltage_label, wavelen_label, current_meter)
        self.play(light_anim)
        self.wait(5)
        self.play(meter_anim_pos)
        self.wait(4)
        self.play(LaggedStart(experiment.voltage.animate(run_time = 5, rate_func = linear).set_value(final_voltage), 
                              current_meter.current.animate(run_time = final_intensity * 7, rate_func = smooth).set_value(0),
                              lag_ratio= 0.4)
                  )
        self.wait(30)
        

class ExperimentDataTable(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = (0.05,0.05,0.05,1)     #type: ignore
        #self.camera.background_rgba = (1,1,1,0)     #type: ignore
        self.camera.frame.scale(1.8)
        #timing and other Variables
        graph_writing_time_real = 2
        graph_writing_time_theory = 1
        
        """ //TODO MAKE BETTER

        Order of metal operations idicies 
                    "Cs": 2.1,
                    "Al": 4.08,
                    "U": 3.6,
                    "Zn" : 4.3,
                    "Au": 5.1,
        """          


        experiment_data = ExperimentData()
        self.add(experiment_data)
        
        cs_anim_real = Write(experiment_data.Cs_graph_real, run_time = graph_writing_time_real, rate_func = smooth)
        al_anim_real = Write(experiment_data.Al_graph_real, run_time = graph_writing_time_real, rate_func = smooth)
        u_anim_real  = Write(experiment_data.U_graph_real , run_time = graph_writing_time_real, rate_func = smooth)
        zn_anim_real = Write(experiment_data.Zn_graph_real, run_time = graph_writing_time_real, rate_func = smooth)
        au_anim_real = Write(experiment_data.Au_graph_real, run_time = graph_writing_time_real, rate_func = smooth)

        cs_anim_theory = Write(experiment_data.Cs_graph_theory, run_time = graph_writing_time_theory, rate_func = smooth)
        al_anim_theory = Write(experiment_data.Al_graph_theory, run_time = graph_writing_time_theory, rate_func = smooth)
        u_anim_theory  = Write(experiment_data.U_graph_theory , run_time = graph_writing_time_theory, rate_func = smooth)
        zn_anim_theory = Write(experiment_data.Zn_graph_theory, run_time = graph_writing_time_theory, rate_func = smooth)
        au_anim_theory = Write(experiment_data.Au_graph_theory, run_time = graph_writing_time_theory, rate_func = smooth)
        
        self.play(LaggedStart(cs_anim_theory, 
                              FadeIn(experiment_data.metal_labels[0].scale(1.5)),
                              lag_ratio= 0.8
                              ))
        
        self.wait(2)
        
        self.play(LaggedStart(
            al_anim_theory, u_anim_theory, zn_anim_theory, au_anim_theory,
            cs_anim_real,al_anim_real,u_anim_real,zn_anim_real,au_anim_real,FadeIn(experiment_data.metal_labels[1:5]),
            lag_ratio= 0.5
        ))
        

class PhotonIntro(Scene):
     def construct(self) -> None:
        self.camera.background_rgba = (0.05,0.05,0.05,1)     #type: ignore
        self.camera.frame.scale(0.5)
        axes = Axes(x_range=(-5,5,1),
                  y_range= (-5,5,1),
                 )
        def t_updater(mob, dt):
             mob.time += dt
        photon = WavePacket2D(axes= axes, wave_len= 400).add_updater(t_updater)
        self.add(photon)
        self.wait(60)


class Electron(VGroup):
    def __init__(
            self, 
            *vmobjects, 
            **kwargs
            ):
        super().__init__(**kwargs)
        main = Circle(stroke_color= BLUE, fill_color = BLUE, fill_opacity = .5).scale(0.4)
        electron_text = Text(text=" e⁻", font_size= 40)

        self.add(main, electron_text)

        self.time = 0
        self.velocity = np.array([0.0,0.0,0.0])
        self.acceleration = np.array([0.0,0.0,0.0])
        
class Electroscope(Mobject):
    def __init__(
            self,
            initial_charge = 0.0,
            scale = 1.0,      
            **kwargs,
    ):
        super().__init__(**kwargs)
        self.charge = ValueTracker(initial_charge)
        self.scale_factor = scale

        def needle_upd(needle):
            needle.rotate(-(abs(self.charge.get_value())/2 * PI), Z_AXIS, needle.get_center())

        radius = self.scale_factor
        self.body = ImageMobject(filename="C:\\Users\\Priyam\\ManimGL\\manim\\Own\\GFS\\assets\\Electroscope_Body.png").scale(radius)#type: ignore
        self.needle= always_redraw(lambda: Line(start=self.body.get_center() + LEFT * 0.05 * radius + DOWN * 0.15 * radius + DOWN*radius* 0.95, #type:ignore
                                                end= self.body.get_center() + LEFT * 0.05 * radius + DOWN * 0.15 * radius + UP * radius * 0.95, #type:ignore
                                                stroke_width = 6 * radius 
                                                ).set_color(WHITE))
        joint = Dot(self.needle.get_center()).scale(0.7).set_color(WHITE)
        self.needle.add_updater(needle_upd)
        self.add( self.needle,self.body, joint)

class Amperemeter(Mobject):
    def __init__(
            self,
            init_current=.0,
            scale = 1.0,
            **kwargs
    ):
        self.scale_fac = scale
        self.current = ValueTracker(init_current)
        
        super().__init__(**kwargs)

        
        
        def needle_upd(needle):
            needle.rotate(-(np.clip(-0.9, 0.9, (self.current.get_value() * 1.8 - 0.9))/2 * PI), Z_AXIS, needle.get_start())

        self.body = ImageMobject(filename="C:\\Users\\Priyam\\ManimGL\\manim\\Own\\GFS\\assets\\current_meter_body.png").scale(self.scale_fac)#type: ignore
        self.needle= always_redraw(lambda: Line(start=self.body.get_center()   + DOWN*self.scale_fac* 1.5, #type:ignore
                                                end= self.body.get_center()  + UP * self.scale_fac * 0.95, #type:ignore
                                                stroke_width = 10 * self.scale_fac 
                                                ).set_color(BLACK))
        
        self.needle.add_updater(needle_upd)
        self.joint = Dot(self.needle.get_start()).scale(1.7 * self.scale_fac).set_color(BLACK)
        self.add(self.body, self.needle, self.joint)
        self.add(Text("A").set_color(GREY_D).scale(4 * self.scale_fac).move_to(self.body))