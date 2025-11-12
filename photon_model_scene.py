from manimlib import *#type: ignore
import numpy as np
import sys, os
sys.path.append(os.path.dirname("C:\\Users\\Priyam\\ManimGL\\manim\\Own\\GFS"))

from GFS.photoelectric_effect import *


class TestScene(Scene):
    def construct(self) -> None:
        self.camera.background_rgba = (0.05,0.05,0.05,0)     #type: ignore
        wave_len = 500

        e_collector_top   = Collector(height= 1.2).shift(2* LEFT + 1 * UP)#type: ignore
        e_collector_bottom = Collector(height= 1.2).shift(2* LEFT + 1 * DOWN)#type: ignore

        electron_collectors = VGroup(e_collector_top, e_collector_bottom)

        collector_and_emitter = PhotonCollectorAndEmitter(collectors= electron_collectors,incoming_light_wave_len= wave_len).shift(2*RIGHT)#type:ignore
        light_source = PhotonEmitter(collectors=[collector_and_emitter],
                                     width= 0.5,
                                     height=0.5,
                                     emitting_area= 0.25,
                                     particle_lifetime= 6,
                                     wave_len= wave_len

                                     ).shift(4*LEFT)#type:ignore

        self.add(light_source, collector_and_emitter, electron_collectors)

class PhotonExperiment(Scene):
    def construct(self) -> None:
        
        experiment = PhotoElectricExperimentPhoton()
        self.add(experiment)


class PhotoElectricExperimentPhoton(VGroup):
    def __init__(
            self,
            incoming_light_wave_len = 350.0,
            incoming_light_intensity = 1.0,
            emitter_metal = "Cs",  
            *vmobjects, 
            **kwargs):
        super().__init__(*vmobjects, **kwargs)
        self.incoming_light_wave_len = incoming_light_wave_len
        self.incoming_light_intensity = incoming_light_intensity
        self.metals = { # work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
                    "Cs": 2.1,
                    "Al": 4.08,
                    "U": 3.6,
                    "Zn" : 4.3,
                    "Au": 5.1,
                    }
        self.emitter_metal = emitter_metal

        e_collector_top   = Collector(height= 1.2).shift(2* LEFT + 1 * UP)#type: ignore
        e_collector_bottom = Collector(height= 1.2).shift(2* LEFT + 1 * DOWN)#type: ignore

        electron_collectors = VGroup(e_collector_top, e_collector_bottom)

        collector_and_emitter = PhotonCollectorAndEmitter(collectors= electron_collectors,incoming_light_wave_len= self.incoming_light_wave_len).shift(2*RIGHT)#type:ignore
        light_source = PhotonEmitter(collectors=[collector_and_emitter],
                                     width= 0.5,
                                     height=0.5,
                                     emitting_area= 0.25,
                                     particle_lifetime= 6,
                                     wave_len= self.incoming_light_wave_len

                                     ).shift(4*LEFT)#type:ignore

        self.add(light_source, collector_and_emitter, electron_collectors)

class PhotonEmitter(VMobject):
    def __init__(
            self,
            collectors,
            width= .2,
            height= 2.0,
            particle_x_boundary = (-10, 10),
            particle_y_boundary = (-10, 10),
            emission_intensity = 0.5,
            wave_len = 400.0,
            particle_speed_range = (1.0, 1.0),
            emitting_area = 0.25,
            particle_lifetime = 5.0,
            velocity_dir_range = np.array([np.array([5.0, 0.3, 0]), np.array([5.0, -0.3, 0])]),
            emission_side = RIGHT,
            emitter_color = WHITE,
            emitt = True,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.width = width
        self.height = height
        self.emission_intensity = emission_intensity
        self.particle_speed_range = particle_speed_range
        self.particle_lifetime = particle_lifetime
        self.emission_side = emission_side
        self.emitt = emitt
        self.velocity_dir_range = np.array([np.array([velocity_dir_range[0][0] * self.emission_side[0], velocity_dir_range[0][1], 0]), np.array([velocity_dir_range[1][0] * self.emission_side[0], velocity_dir_range[1][1], 0])])
        self.particle_boundary = (particle_x_boundary, particle_y_boundary)
        self.collectors = collectors
        self.wave_len = wave_len

        emitter_body = RoundedRectangle(width= self.width, height= self.height, corner_radius= 0.1 * self.width, fill_color = emitter_color, fill_opacity = 1, stroke_color = emitter_color, **kwargs).set_opacity(0)
        self.emissive_range = (emitter_body.get_y()+ 1/2 * emitting_area, emitter_body.get_y()- 1/2 * emitting_area)

        
        def mov_upd(p, dt):
            p.shift(p.velocity * dt)
            p.time += dt
            

        def destructor(p_group):
            for p in p_group:
                p_x = p.get_x()
                p_y = p.get_y()
                if p.time > self.particle_lifetime or p_x < self.particle_boundary[0][0] or p_x > self.particle_boundary[0][1] or p_y < self.particle_boundary[1][0] or p_y > self.particle_boundary[1][1]:
                    p_group.remove(p)

        def coll_upd(p_group):
            for p in p_group:
                for collector  in self.collectors:
                    if collector.isInside(p.get_center()):
                        p_group.remove(p)

        self.photon_group = VGroup().add_updater(destructor).add_updater(coll_upd)

        def particle_spawner(emitter, dt):
            if self.emitt:
                for _ in range(int(self.emission_intensity* 100 * dt)) :
                    p_x = emitter_body.get_center()[0]
                    p_y = random.uniform(self.emissive_range[0], self.emissive_range[1])
                    p_v = normalize(np.array([random.uniform(self.velocity_dir_range[0][0], self.velocity_dir_range[1][0]), random.uniform(velocity_dir_range[0][1],velocity_dir_range[1][1]), 0])) * random.uniform(self.particle_speed_range[0], self.particle_speed_range[1])#type: ignore                       
                    emitter.photon_group.add(EmittedPhoton(point=([p_x, p_y, 0]), fill_color= self.get_light_color(self.wave_len), #type:ignore
                                                                v= p_v).add_updater(mov_upd))
                    
        self.add_updater(particle_spawner)

        main = Circle(stroke_color= WHITE, fill_color = self.get_light_color(self.wave_len), fill_opacity = self.emission_intensity*0.5, stroke_opacity = 1 ).scale(0.4)
        cross= Cross(main, stroke_color = WHITE).scale(0.7)

        self.add( self.photon_group, emitter_body, main, cross)

    def get_light_color(self, wavelength):
        
        """
        Maps a wavelength (arbitrary units, e.g. ~200 nm per unit) to a Manim RGB color.
        Below min_wl → violet
        Above max_wl → red
        Between → smooth gradient.
        Returns a Color created from RGB values in [0, 1].
        """

        # Define the wavelength limits (in your arbitrary units)
        min_wl = 380
        max_wl = 800

        wl = np.clip(wavelength, min_wl, max_wl)

        # Convert to nanometers (roughly)
        wl_nm = np.interp(wl, [min_wl, max_wl], [380, 700])

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

        R = np.clip(R * factor, 0, 1)
        G = np.clip(G * factor, 0, 1)
        B = np.clip(B * factor, 0, 1)

        # Return as Manim Color
        return Color(rgb=(R, G, B))
        
class PhotonCollectorAndEmitter(VMobject):
    def __init__(
            self,
            collectors= VGroup(), # provide a list with all collectors
            width= .2,
            height= 3.0,
            incoming_light_wave_len = 400.0,
            metal = "Cs",
            particle_x_boundary = (-10, 2),
            particle_y_boundary = (-10, 10),
            emission_intensity = 0.5,
            particle_lifetime = 15.0,
            velocity_dir_range = np.array([np.array([5.0, 1.0, 0]), np.array([5.0, -1.0, 0])]),
            emission_side = LEFT,
            emitter_color = WHITE,
            emitt = True,
            **kwargs
    ):
        self.width = width
        self.height = height
        self.emitter_color = emitter_color
        self.metal = { # work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
                    "Cs": 2.1,
                    "Al": 4.08,
                    "U": 3.6,
                    "Zn" : 4.3,
                    "Au": 5.1,
                    }
        self.photon_collision_counter = 0
        self.emission_intensity = emission_intensity
        self.collectors = collectors
        self.particle_lifetime = particle_lifetime
        self.emission_side = emission_side
        self.emitt = emitt
        self.velocity_dir_range = np.array([np.array([velocity_dir_range[0][0] * self.emission_side[0], velocity_dir_range[0][1], 0]), np.array([velocity_dir_range[1][0] * self.emission_side[0], velocity_dir_range[1][1], 0])])
        self.particle_boundary = (particle_x_boundary, particle_y_boundary)
        self.incoming_light_wave_len = incoming_light_wave_len

        self.particle_speed_range = (self.wave_len_to_max_speed(wave_len=self.incoming_light_wave_len, work_func= self.metal[metal])/200000 -1, self.wave_len_to_max_speed(wave_len=self.incoming_light_wave_len, work_func= self.metal[metal])/200000)
        
        if self.get_stopping_voltage(self.incoming_light_wave_len, self.metal[metal]) < 0:
            self.emission_intensity = 0
       
        super().__init__(**kwargs)
        
        self.body = RoundedRectangle(width= self.width, 
                                    height= self.height, 
                                    fill_color= self.emitter_color, 
                                    fill_opacity = 1,
                                    stroke_color = self.emitter_color,
                                    stroke_opacity = 0,
                                    corner_radius= 0.1 * self.width, 
                                    **kwargs)
        
        
        
        def destructor(e_group):
            for e in e_group:
                e_x = e.get_x()
                e_y = e.get_y()
                if e.time > self.particle_lifetime or e_x < self.particle_boundary[0][0] or e_x > self.particle_boundary[0][1] or e_y < self.particle_boundary[1][0] or e_y > self.particle_boundary[1][1]:
                    e_group.remove(e)
        
        def t_upd(e_group, dt):
            for e in e_group:
                e.time += dt
        

        def coll_upd(e_group):
            for e in e_group:
                if len(self.collectors) > 0:
                    for collector in self.collectors:
                        if collector.isInside(e.get_center()): # type: ignore
                            e_group.remove(e)
        
        self.electron_group = VGroup().add_updater(t_upd).add_updater(destructor).add_updater(coll_upd)
    
        self.add(self.electron_group ,self.body)
    
    def mov_upd(self, e, dt):
            e.shift(e.velocity * dt)
            e.time += dt
    
    def isInside(self, point):
        out = False

        if point[0] < self.body.get_edge_center(RIGHT)[0] and point[0] > self.body.get_edge_center(LEFT)[0] and point[1] < self.body.get_edge_center(UP)[1] and point[1] > self.body.get_edge_center(DOWN)[1]:
            out = True
            self.photon_collision_counter += 1
            p_v = normalize(np.array([random.uniform(self.velocity_dir_range[0][0], self.velocity_dir_range[1][0]), random.uniform(self.velocity_dir_range[0][1],self.velocity_dir_range[1][1]), 0])) * random.uniform(self.particle_speed_range[0], self.particle_speed_range[1])#type: ignore
            self.electron_group.add(EmittedElectron(point= point, v= p_v).add_updater(self.mov_upd))
        return out
    
    def wave_len_to_max_speed(self, wave_len, work_func):
        #wave_len to freq
        c= 299792458
        true_wave_len = wave_len / 1000000000
        freq = c/true_wave_len

        electron_mass = 9.109e-31
        electron_charge=1.602e-19

        
        
        speed = np.sqrt((2*electron_charge*self.get_stopping_voltage(freq, work_func))/electron_mass)
        
        return speed

    def get_stopping_voltage(self, frequency, work_func):
            m = 4.1357e-15
            voltage = m*frequency - work_func
            return voltage
    
    
class EmittedPhoton(Dot):
    def __init__(self, 
                 point = ORIGIN, 
                 radius = .01, 
                 stroke_color = WHITE, 
                 stroke_width = 0, 
                 fill_opacity = 1, 
                 fill_color = YELLOW_A,
                 v = np.array([0,0,0]),
                 **kwargs
                 ):
        super().__init__(point, 
                         radius, 
                         stroke_color, 
                         stroke_width, 
                         fill_opacity, 
                         fill_color, 
                         **kwargs)
        
        self.velocity = v
        self.time = 0