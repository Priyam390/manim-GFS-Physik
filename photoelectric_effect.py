
from manimlib import * #type: ignore
import numpy as np
import random 

#work Functions in eV for different metals
# work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
WF_CS = 1.94 # origninal: 2.1
WF_AL = 4.08# origninal: 4.08
WF_U  = 3.6# origninal: 3.6
WF_ZN = 4.27# origninal: 4.3
WF_AU = 5.1# origninal: 5.1



class PhotoElectricExperimentBasic(VGroup):
    def __init__(
            self,
            incoming_light_wave_len = 250.0,
            incoming_light_intensity = 1.0,
            always_redraw_e_field = False,
            emitter_metal = "Cs",
            max_voltage=3.0,
            performant_eField = True,
            simulation_speed = 1.0,
            electron_lifetime = 7,
            electron_x_boundary=(-5,2),
            electron_y_boundary=(-3,3),
            *vmobjects,
            **kwargs
            ):
        
        super().__init__(*vmobjects, **kwargs)
        
        self.incoming_light_wave_len = incoming_light_wave_len
        self.incoming_light_intensity = incoming_light_intensity
        self.emitter_metal = emitter_metal
        self.emission_intensity = incoming_light_intensity
        self.sim_speed = simulation_speed
        self.electron_lifetime = electron_lifetime
        self.electron_x_boundary = electron_x_boundary
        self.electron_y_boundary = electron_y_boundary

        self.metals = { # work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
                    "Cs": WF_CS,
                    "Al": WF_AL,
                    "U": WF_U,
                    "Zn" : WF_ZN,
                    "Au": WF_AU,
                    }

        self.axes = NumberPlane()

        self.voltage = ValueTracker(0)

        def voltage_upd(e_field):
            e_field.voltage  = -self.voltage.get_value()

        self.lightsource = LightCone(wave_length= self.incoming_light_wave_len, intensity= self.incoming_light_intensity)
        self.add(self.lightsource)

        if always_redraw_e_field:
            electric_field_calc =always_redraw(lambda: PlateCapacitorEField(edges=[-1.9 , 2, -1.5, 1.5], coordinate_system= self.axes, voltage= self.voltage.get_value()).add_updater(voltage_upd))#type: ignore
            self.add(electric_field_calc)#type: ignore
        elif not always_redraw_e_field: 
            electric_field_calc =PlateCapacitorEField(edges=[-1.9 , 2, -1.5, 1.5], coordinate_system= self.axes, voltage= self.voltage.get_value()).add_updater(voltage_upd)#type: ignore
            self.add(electric_field_calc)#type: ignore
        
        if performant_eField:
            self.show_e_field = ShowEField(max_voltage=max_voltage).add_updater(voltage_upd)
            self.add(self.show_e_field)

        if self.wave_len_to_max_speed(incoming_light_wave_len, self.metals[self.emitter_metal]) <= 0:
            self.emission_intensity = 0

        max_speed = (self.wave_len_to_max_speed(incoming_light_wave_len, self.metals[self.emitter_metal])/205000)

        particle_speed_calc = (max_speed *0.6, max_speed)

        self.collector_top = Collector(height= 1.2).shift(2* LEFT + 1 * UP)#type: ignore
        self.collector_bottom = Collector(height= 1.2).shift(2* LEFT + 1 * DOWN)#type: ignore
        
        self.emitter = Emitter(emission_intensity= self.emission_intensity,
                          height= 3, 
                          velocity_dir_range= np.array([np.array([1.5, 1.0, 0]), np.array([1.5, -1.0, 0])]),
                          particle_speed_range= particle_speed_calc, 
                          emitting_area= 1,
                          particle_lifetime= self.electron_lifetime,
                          particle_x_boundary=self.electron_x_boundary,
                          particle_y_boundary=self.electron_y_boundary,
                          emission_side= LEFT,
                          emitter_color= GREY_D,
                          electric_field= True,
                          e_field= electric_field_calc,
                          simulation_speed= self.sim_speed,
                          collectors= [self.collector_top, self.collector_bottom]).shift(2*RIGHT)#type: ignore
        self.emitter.emitt = True

        self.add(self.collector_top, self.collector_bottom, self.emitter)
        self.total_collisions = 0
        def coll_upd(mob):
            mob.total_collisions = self.collector_top.collision_counter + self.collector_bottom.collision_counter

        self.add_updater(coll_upd)
    
        self.electron_speed_real = self.wave_len_to_max_speed(incoming_light_wave_len, self.metals["Cs"])
        
        

    def wave_len_to_max_speed(self, wave_len, work_func):
        #wave_len to freq
        c= 299792458
        true_wave_len = wave_len / 1000000000
        freq = c/true_wave_len

        electron_mass = 9.109e-31
        electron_charge=1.602e-19

        def get_stopping_voltage(frequency, work_func):
            m = 4.1357e-15
            voltage = m*frequency - work_func
            return voltage
        
        speed = np.sqrt((2*electron_charge*get_stopping_voltage(freq, work_func))/electron_mass)
        
        return speed

class LightCone(VMobject):
    def __init__(
            self,
            wave_length,
            intensity = 0.3,
            auto_create = True,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.intensity = intensity
        self.wave_len = wave_length
        self.cone = always_redraw(lambda :Polygon([-5.65,-0.25,0],[-5.55,0,0],[-5.65,0.25,0],[2,0.5,0],[2, -0.5, 0], stroke_opacity = 0, fill_color = self.get_light_color(self.wave_len), fill_opacity = self.intensity * 0.5))# type: ignore
 
        main = Circle(stroke_color= WHITE, fill_color = self.get_light_color(self.wave_len), fill_opacity = intensity*0.5, stroke_opacity = 1 ).scale(0.4)
        cross= Cross(main, stroke_color = WHITE).scale(0.7)
        
        self.source = VGroup(main, cross).shift(6 * LEFT)#type: ignore

        if auto_create:
            self.add(self.cone, self.source)#type: ignore

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

        if wavelength < 380 or wavelength > 700:
            return Color(rgb=(0.2, 0.2, 0.2))


        # Convert to nanometers (roughly)
        wl_nm = np.interp(wavelength, [min_wl, max_wl], [380, 700])

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

class Emitter(VMobject):
    def __init__(
            self,
            collectors,
            width= .2,
            height= 2.0,
            particle_x_boundary = (-10, 10),
            particle_y_boundary = (-10, 10),
            emission_intensity = 1.0,
            particle_speed_range = (1.0, 1.0),
            emitting_area = 0.75,
            particle_lifetime = 5.0,
            velocity_dir_range = np.array([np.array([5.0, 1.0, 0]), np.array([5.0, -1.0, 0])]),
            emission_side = RIGHT,
            emitter_color = WHITE,
            electric_field = False,
            e_field = None,
            emitt = False,
            simulation_speed = 1.0,
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
        self.sim_speed = simulation_speed

        self.emitter_body = RoundedRectangle(width= self.width, height= self.height, corner_radius= 0.1 * self.width, fill_color = emitter_color, fill_opacity = 1, stroke_color = emitter_color, **kwargs)
        self.emissive_range = (self.emitter_body.get_y()+ 1/2 * emitting_area, self.emitter_body.get_y()- 1/2 * emitting_area)

        if not electric_field:
            def mov_upd(e, dt):
                e.shift(e.velocity * dt * self.sim_speed)
                e.time += dt * self.sim_speed
            

            def destructor(e_group):
                for e in e_group:
                    e_x = e.get_x()
                    e_y = e.get_y()
                    if e.time > self.particle_lifetime or e_x < self.particle_boundary[0][0] or e_x > self.particle_boundary[0][1] or e_y < self.particle_boundary[1][0] or e_y > self.particle_boundary[1][1]:
                        e_group.remove(e)

            def coll_upd(e_group):
                for e in e_group:
                    for collector  in self.collectors:
                        if collector.isInside(e.get_center()):
                            e_group.remove(e)

            self.electron_group = VGroup().add_updater(destructor).add_updater(coll_upd)

            def particle_spawner(emitter, dt):
                if self.emitt:
                    for _ in range(int(self.emission_intensity* 100 * dt * self.sim_speed)) :
                        
                        p_x = self.emitter_body.get_center()[0]

                        p_y = random.uniform(self.emissive_range[0], self.emissive_range[1])

                        p_v = normalize(np.array([random.uniform(emitter.velocity_dir_range[0][0], emitter.velocity_dir_range[1][0]), random.uniform(velocity_dir_range[0][1],velocity_dir_range[1][1]), 0])) * random.uniform(self.particle_speed_range[0], self.particle_speed_range[1])#type: ignore
                        
                        emitter.electron_group.add(EmittedElectron(point=([p_x, p_y, 0]), #type:ignore
                                                                v= p_v).add_updater(mov_upd))
                    
            self.add_updater(particle_spawner)

            self.add( self.electron_group, self.emitter_body)
        
        elif electric_field:
            def mov_upd(e, dt):
                e.shift(e.velocity * dt * self.sim_speed)
                e.time += dt * self.sim_speed

            def velocity_upd(e, dt):
                
                force = 1.602e-19 * (e_field.vect_func([e.get_center()])[0]/0.02)#type: ignore
                electron_mass = electron_mass = 9.109e-31
                
                e.velocity += ((force/electron_mass) / 8000000000000) * dt * self.sim_speed #type: ignore


            def destructor(e_group):
                for e in e_group:
                    e_x = e.get_x()
                    e_y = e.get_y()
                    if e.time > self.particle_lifetime or e_x < self.particle_boundary[0][0] or e_x > self.particle_boundary[0][1] or e_y < self.particle_boundary[1][0] or e_y > self.particle_boundary[1][1]:
                        e_group.remove(e)
            
            def coll_upd(e_group):
                for e in e_group:
                    for collector  in self.collectors:
                        if collector.isInside(e.get_center()):
                            e_group.remove(e)
            
            self.electron_group = VGroup().add_updater(destructor).add_updater(coll_upd)

            def particle_spawner(emitter, dt):
                if self.emitt:        
                    for _ in range(int(self.emission_intensity* 100 * dt * self.sim_speed)) :
                        
                        p_x = self.emitter_body.get_edge_center(self.emission_side)[0]

                        p_y = random.uniform(self.emissive_range[0], self.emissive_range[1])

                        p_v = normalize(np.array([random.uniform(emitter.velocity_dir_range[0][0], emitter.velocity_dir_range[1][0]), random.uniform(velocity_dir_range[0][1],velocity_dir_range[1][1]), 0])) * random.uniform(self.particle_speed_range[0], self.particle_speed_range[1])#type: ignore
                        
                        emitter.electron_group.add(EmittedElectron(point=([p_x, p_y, 0]), #type:ignore
                                                                v= p_v).add_updater(mov_upd).add_updater(velocity_upd))
                    
            self.add_updater(particle_spawner)

            self.add(self.electron_group, self.emitter_body)
               
class Collector(VMobject):
    def __init__(
            self,
            width= .2,
            height= 2.0,
            collector_color = WHITE,
            
            **kwargs
    ):
        self.width = width
        self.height = height
        self.collector_color = collector_color
        self.collision_counter = 0
       
        super().__init__(**kwargs)
        
        self.body = RoundedRectangle(width= self.width, 
                                    height= self.height, 
                                    fill_color= self.collector_color, 
                                    fill_opacity = 1,
                                    stroke_color = self.collector_color,
                                    stroke_opacity = 0,
                                    corner_radius= 0.1 * self.width, 
                                    **kwargs)
        
        def destructor(m_group):
            for m in m_group:
                if m.time >= m.lifespan:
                    m_group.remove(m)
        
        def t_upd(m_group, dt):
            for m in m_group:
                m.time += dt
        self.hitmarker_group = VGroup().add_updater(t_upd).add_updater(destructor)
        

        self.add(self.body,self.hitmarker_group)
       
    
    def isInside(self, point):
        out = False

        if point[0] < self.body.get_edge_center(RIGHT)[0] and point[0] > self.body.get_edge_center(LEFT)[0] and point[1] < self.body.get_edge_center(UP)[1] and point[1] > self.body.get_edge_center(DOWN)[1]:
            out = True
            self.collision_counter += 1
            self.hitmarker_group.add(Hitmarker(point= point))
        return out

class PlateCapacitorEField(VectorField):
    def __init__(
            self, 
            edges, #in shape [left_edge_x, right_edge_x, bottom_edge_y, top_edge_y] all edges are either y or x coordinates
            coordinate_system,
            voltage = 1.0,
            voltage_dir = LEFT, 
            density= 2, 
            
            color = GREEN, 
            color_map_name= "3b1b_colormap", 
            color_map = None, 
            stroke_opacity= 0.35, 
            stroke_width= 3, 
            tip_width_ratio= 4.0, 
            tip_len_to_width= 0.01, 
            max_vect_len = None, 
            max_vect_len_to_step_size= 0.8, 
            flat_stroke = False, 
            norm_to_opacity_func=None, 
            **kwargs
            ):
        
        self.voltage = voltage
        self.boundary = edges
        self.voltage_dir = voltage_dir
        self.axes = coordinate_system.set_opacity(0)
        
        
    

        super().__init__(func=self.vect_func, #type: ignore
                         coordinate_system=self.axes, 
                         density=density, 
                         color=color, 
                         color_map_name=color_map_name, 
                         color_map=color_map, 
                         stroke_opacity=stroke_opacity, 
                         stroke_width=stroke_width, 
                         tip_width_ratio=tip_width_ratio, 
                         tip_len_to_width=tip_len_to_width, 
                         max_vect_len=max_vect_len, 
                         max_vect_len_to_step_size=max_vect_len_to_step_size, 
                         flat_stroke=flat_stroke, 
                         norm_to_opacity_func=norm_to_opacity_func, 
                         **kwargs)     

        self.add(self.axes)

    def vect_func(self, coords):
        out =[]
        for coord in coords:
            if coord[0] >= self.boundary[0] and coord[0] <= self.boundary[1] and coord[1] >= self.boundary[2] and coord[1] <= self.boundary[3]:
                out.append(self.voltage_dir * self.voltage)#type: ignore
            else:
                out.append(np.array([0,0,0]))
        return np.array(out)
    
class EmittedElectron(Dot):
    def __init__(self, 
                 point = ORIGIN, 
                 radius = .015, 
                 stroke_color = WHITE, 
                 stroke_width = 0, 
                 fill_opacity = 1, 
                 fill_color = BLUE_B,
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
        
class Hitmarker(Dot):
    def __init__(self, 
                 point, 
                 radius = 0.05, 
                 stroke_width = 0,
                 stroke_color = RED_E, 
                 fill_opacity = 0.9, 
                 fill_color = RED_E,
                 lifespan =0.5, 
                 **kwargs):
        
        self.time = 0
        self.lifespan = lifespan
        self.opacity = fill_opacity
        
        super().__init__(point, 
                         radius, 
                         stroke_color, 
                         stroke_width, 
                         self.opacity, 
                         fill_color, 
                         **kwargs)
        
class ExperimentData(VGroup):
    def __init__(
            self, 
            *vmobjects,
            line_width = 2,
            real_opacity = 1,
            theory_opacity = 0.5,
            real_x_range = (-1, 19.5),
            theroy_x_range = (-1, 19.5),
            num_dashes_theory = 150,
            appropriate_color = WHITE,
            **kwargs,
            ):
        super().__init__(*vmobjects, **kwargs)
        self.appropriate_color = appropriate_color

        slope = 4.135667707248e-1 

        self.a = Axes(x_range=(-1,20,1),
                    y_range=(-6,8,1),
                   
                    x_axis_config={
                        "include_numbers" :True,
                        "numbers_to_exclude" : [0,20],
                        
                    },
                        ).set_color(self.appropriate_color)
        self.numberplane = NumberPlane(x_range=(-1,20,1),
                                  y_range=(-6,8,1),
                                  background_line_style={
                                      "stroke_color" : GREY_D,
                                      "stroke_opacity" : 0.5,
                                      "stroke_width": 0.5
                                  }).shift(0.1 * RIGHT) #type: ignore
        self.add_to_back(self.numberplane)
        
        self.metals = { # work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
                    "Cs": WF_CS,
                    "Al": WF_AL,
                    "U": WF_U,
                    "Zn" : WF_ZN,
                    "Au": WF_AU,
                    }
        
        y_axis_nums = VGroup(
            DecimalNumber(number = -5+ i, num_decimal_places= 0, font_size= 25).move_to(self.a.c2p(0,-5 + i)).shift(0.3 * LEFT)# type: ignore
            for i in range(13)
        )

        y_axis_nums[5].shift(0.3 * DOWN)    #type: ignore    

        axis_units = VGroup(
            
            Text(text= "U in V", font_size= 35).move_to(self.a.c2p(0,8)).shift(1.2 * RIGHT), # type: ignore
            Text(text= "f in 10  Hz", font_size= 35).move_to(self.a.c2p(19,-1)), # type: ignore
            Text(text= "14", font_size= 20).move_to(self.a.c2p(19.53,-0.8)), # type: ignore
        )
        
        self.metal_labels = VGroup(
            Text(text= "Cs", font_size= 35).move_to(self.a.c2p(20, slope*20 - self.metals["Cs"])), # type: ignore
            Text(text= "Al", font_size= 35).move_to(self.a.c2p(20, slope*20 - self.metals["Al"])), # type: ignore
            Text(text= "U", font_size= 35).move_to(self.a.c2p(20, slope*20 - self.metals["U"])), # type: ignore
            Text(text= "Zn", font_size= 35).move_to(self.a.c2p(20, slope*20 - self.metals["Zn"] - 0.2)), # type: ignore
            Text(text= "Au", font_size= 35).move_to(self.a.c2p(20, slope*20 - self.metals["Au"])), # type: ignore
        )
        for obj in self.metal_labels:
            obj.set_color(self.appropriate_color)
        for obj in axis_units:
            obj.set_color(self.appropriate_color)
        
        self.add(axis_units, y_axis_nums)

 
        self.line_width = line_width
        self.real_opacity = real_opacity
        self.theory_opacity = theory_opacity
        self.real_x_range = real_x_range
        self.theroy_x_range = theroy_x_range
        self.num_dashes_theory = num_dashes_theory

        self.Cs_graph_real = self.a.get_graph(lambda x: np.clip(slope*x - self.metals["Cs"], 0, 10), stroke_color = self.appropriate_color, stroke_opacity = self.real_opacity, x_range= (-1,self.real_x_range[1]), stroke_width = self.line_width)
        self.Al_graph_real = self.a.get_graph(lambda x: np.clip(slope*x - self.metals["Al"], 0, 10),stroke_color = PURPLE, stroke_opacity = self.real_opacity, x_range= (-1,self.real_x_range[1]), stroke_width = self.line_width)
        self.U_graph_real  = self.a.get_graph(lambda x: np.clip(slope*x - self.metals["U"], 0, 10), stroke_color = GREEN, stroke_opacity = self.real_opacity, x_range= (-1,self.real_x_range[1]), stroke_width = self.line_width)
        self.Zn_graph_real = self.a.get_graph(lambda x: np.clip(slope*x - self.metals["Zn"], 0, 10), stroke_color = BLUE, stroke_opacity = self.real_opacity, x_range= (-1,self.real_x_range[1]), stroke_width = self.line_width)
        self.Au_graph_real = self.a.get_graph(lambda x: np.clip(slope*x - self.metals["Au"], 0, 10), stroke_color = GOLD, stroke_opacity = self.real_opacity, x_range= (-1,self.real_x_range[1]), stroke_width = self.line_width)

        self.Cs_graph_theory = DashedVMobject(self.a.get_graph(lambda x: slope*x - self.metals["Cs"], stroke_color = self.appropriate_color, stroke_opacity = self.theory_opacity,  x_range= self.theroy_x_range, stroke_width = self.line_width), num_dashes= self.num_dashes_theory)
        self.Al_graph_theory = DashedVMobject(self.a.get_graph(lambda x: slope*x - self.metals["Al"],stroke_color = PURPLE, stroke_opacity = self.theory_opacity,  x_range= self.theroy_x_range, stroke_width = self.line_width), num_dashes= self.num_dashes_theory)
        self.U_graph_theory  = DashedVMobject(self.a.get_graph(lambda x: slope*x - self.metals["U"], stroke_color = GREEN, stroke_opacity = self.theory_opacity,  x_range= self.theroy_x_range, stroke_width = self.line_width), num_dashes= self.num_dashes_theory)
        self.Zn_graph_theory = DashedVMobject(self.a.get_graph(lambda x: slope*x - self.metals["Zn"], stroke_color = BLUE, stroke_opacity = self.theory_opacity,  x_range= self.theroy_x_range, stroke_width = self.line_width), num_dashes= self.num_dashes_theory)
        self.Au_graph_theory = DashedVMobject(self.a.get_graph(lambda x: slope*x - self.metals["Au"], stroke_color = GOLD, stroke_opacity = self.theory_opacity,  x_range= self.theroy_x_range, stroke_width = self.line_width), num_dashes= self.num_dashes_theory)

        #self.add(self.Cs_graph_theory, self.Al_graph_theory, self.U_graph_theory, self.Zn_graph_theory, self.Au_graph_theory)

        #self.add(self.Cs_graph_real, self.Al_graph_real, self.U_graph_real, self.Zn_graph_real, self.Au_graph_real)
        self.add(self.a)

class PhotoElectricExperimentPhoton(VGroup):
    def __init__(
            self,
            incoming_light_wave_len = 420.0,
            incoming_light_intensity = 0.5,
            emitter_metal = "Cs", 
            voltage = 0, 
            always_redraw_e_field = False,
            max_voltage = 5.0,
            performant_eField = True,
            electron_lifetime = 7,
            simulation_speed = 1,
            *vmobjects, 
            **kwargs):
        super().__init__(*vmobjects, **kwargs)
        self.incoming_light_wave_len = incoming_light_wave_len
        self.incoming_light_intensity = incoming_light_intensity
        self.metals = { # work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
                    "Cs": WF_CS,
                    "Al": WF_AL,
                    "U": WF_U,
                    "Zn" : WF_ZN,
                    "Au": WF_AU,
                    }
        self.emitter_metal = emitter_metal
        self.sim_speed = simulation_speed
        self.voltage = ValueTracker(voltage * self.sim_speed)
        self.axes = NumberPlane()
        

        

        e_collector_top   = Collector(height= 1.2).shift(2* LEFT + 1 * UP)#type: ignore
        e_collector_bottom = Collector(height= 1.2).shift(2* LEFT + 1 * DOWN)#type: ignore

        self.electron_collectors = VGroup(e_collector_top, e_collector_bottom)

        def voltage_upd(e_field):
            e_field.voltage = self.voltage.get_value()
        
        
        if always_redraw_e_field:
            electric_field = always_redraw(lambda : PlateCapacitorEField(edges=[-1.9 , 2, -1.5, 1.5],coordinate_system= self.axes, voltage= self.voltage.get_value(),).add_updater(voltage_upd))#type: ignore
            self.add(electric_field)#type:ignore
        elif not always_redraw_e_field:
            electric_field = PlateCapacitorEField(edges=[-1.9 , 2, -1.5, 1.5],coordinate_system= self.axes, voltage= self.voltage.get_value(),).add_updater(voltage_upd)#type: ignore
            self.add(electric_field)#type:ignore
        if performant_eField:
            self.show_e_field = ShowEField(max_voltage=max_voltage).add_updater(voltage_upd)
            self.add(self.show_e_field)

        self.collector_and_emitter = PhotonCollectorAndEmitter(e_field=True, 
                                                          e_field_obj= electric_field,
                                                          collectors= self.electron_collectors,
                                                          incoming_light_wave_len= self.incoming_light_wave_len,
                                                          particle_lifetime= electron_lifetime,
                                                          simulation_speed= self.sim_speed
                                                          ).shift(2*RIGHT)#type:ignore
        

        self.light_source = PhotonEmitter(
                                     collectors=[self.collector_and_emitter],
                                     width= 0.5,
                                     height=0.5,
                                     emitting_area= 0.1,
                                     particle_lifetime= 3.5,
                                     wave_len= self.incoming_light_wave_len,
                                     particle_speed_range=(2,2),
                                     emission_intensity= self.incoming_light_intensity,
                                     simulation_speed= self.sim_speed,

                                     ).shift(4*LEFT)#type:ignore

        self.add(self.light_source, self.collector_and_emitter, self.electron_collectors)
  

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
            velocity_dir_range = np.array([np.array([5.0, 0.5, 0]), np.array([5.0, -0.5, 0])]),
            emission_side = RIGHT,
            emitter_color = WHITE,
            emitt = True,
            simulation_speed = 1.0,
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
        self.sim_speed = simulation_speed

        main =always_redraw(lambda: Circle(stroke_color= WHITE, fill_color = self.get_light_color(self.wave_len), fill_opacity = self.emission_intensity*0.5, stroke_opacity = 1 ).scale(0.4).shift(4*LEFT))
        cross= always_redraw(lambda: Cross(main, stroke_color = WHITE).scale(0.7))


        self.emitter_body = RoundedRectangle(width= self.width, height= self.height, corner_radius= 0.1 * self.width, fill_color = emitter_color, fill_opacity = 1, stroke_color = emitter_color, **kwargs).set_opacity(0).shift(4*RIGHT)
        self.emissive_range = (self.emitter_body.get_y()+ 1/2 * emitting_area, self.emitter_body.get_y()- 1/2 * emitting_area)
        self.add(self.emitter_body)
        
        self.center_coords = self.get_center()

        def mov_upd(p, dt):
            p.shift(p.velocity * dt * self.sim_speed)
            p.time += dt * self.sim_speed
            

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
        self.time_counter = 0
        self.treshhold = 0.017 * 5

        def particle_spawner(emitter, dt):
            self.time_counter += dt * self.sim_speed
            if self.emitt and not self.treshhold:
                for _ in range(int(self.emission_intensity* 5)) :
                    p_x = random.uniform(self.center_coords[0] - emitting_area/2, self.center_coords[0]+ emitting_area/2)
                    p_y = random.uniform(self.emissive_range[0], self.emissive_range[1])
                    p_v = normalize(np.array([random.uniform(self.velocity_dir_range[0][0], self.velocity_dir_range[1][0]), random.uniform(velocity_dir_range[0][1],velocity_dir_range[1][1]), 0])) * random.uniform(self.particle_speed_range[0], self.particle_speed_range[1])#type: ignore                       
                    emitter.photon_group.add(EmittedPhoton(point=([p_x, p_y, 0]), fill_color= self.get_light_color(self.wave_len), #type:ignore
                                                                v= p_v).add_updater(mov_upd))
                    self.time_counter = 0
            if self.emitt and self.treshhold:
                if self.treshhold - self.time_counter < 0:
                    for _ in range(int(self.emission_intensity* 5)) :
                        p_x = random.uniform(self.center_coords[0] - emitting_area/2, self.center_coords[0]+ emitting_area/2)
                        p_y = random.uniform(self.emissive_range[0], self.emissive_range[1])
                        p_v = normalize(np.array([random.uniform(self.velocity_dir_range[0][0], self.velocity_dir_range[1][0]), random.uniform(velocity_dir_range[0][1],velocity_dir_range[1][1]), 0])) * random.uniform(self.particle_speed_range[0], self.particle_speed_range[1])#type: ignore                       
                        emitter.photon_group.add(EmittedPhoton(point=([p_x, p_y, 0]), fill_color= self.get_light_color(self.wave_len), #type:ignore
                                                                    v= p_v).add_updater(mov_upd))
                    self.time_counter = 0
                    
        self.add_updater(particle_spawner)

        

        self.add( self.photon_group,  main, cross)

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
            e_field = False,
            e_field_obj = None,
            simulation_speed = 1.0,
            **kwargs
    ):
        self.width = width
        self.height = height
        self.emitter_color = emitter_color
        self.metal = { # work functions of different metals in eV from : http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
                    "Cs": WF_CS,
                    "Al": WF_AL,
                    "U": WF_U,
                    "Zn" : WF_ZN,
                    "Au": WF_AU,
                    }
        self.plate_material = metal
        self.photon_collision_counter = 0
        self.emission_intensity = emission_intensity
        self.collectors = collectors
        self.particle_lifetime = particle_lifetime
        self.emission_side = emission_side
        self.emitt = emitt
        self.velocity_dir_range = np.array([np.array([velocity_dir_range[0][0] * self.emission_side[0], velocity_dir_range[0][1], 0]), np.array([velocity_dir_range[1][0] * self.emission_side[0], velocity_dir_range[1][1], 0])])
        self.particle_boundary = (particle_x_boundary, particle_y_boundary)
        self.incoming_light_wave_len = incoming_light_wave_len
        self.is_in_e_field = e_field
        self.e_field = e_field_obj
        self.sim_speed = simulation_speed

        self.particle_speed_range = ((self.wave_len_to_max_speed(wave_len=self.incoming_light_wave_len, work_func= self.metal[self.plate_material])/200000)*0.6, self.wave_len_to_max_speed(wave_len=self.incoming_light_wave_len, work_func= self.metal[self.plate_material])/208000)
        
        if self.get_stopping_voltage(self.incoming_light_wave_len, self.metal[self.plate_material]) < 0:
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
                e.time += dt* self.sim_speed
        

        def coll_upd(e_group):
            for e in e_group:
                if len(self.collectors) > 0:
                    for collector in self.collectors:
                        if collector.isInside(e.get_center()): # type: ignore
                            e_group.remove(e)
        def v_updater(e_group, dt):              
            electron_mass = 9.109e-31   
            for e in e_group:
                force = 1.602e-19 * (self.e_field.vect_func([e.get_center()])[0]/0.02)#type: ignore 
                e.velocity -= ((force/electron_mass) / 8000000000000) * dt * self.sim_speed #type: ignore
        
        self.electron_group = VGroup().add_updater(t_upd).add_updater(destructor).add_updater(coll_upd)

        if self.is_in_e_field:
            self.electron_group.add_updater(v_updater)

        self.add(self.electron_group ,self.body)
    
    def mov_upd(self, e, dt):
            e.shift(e.velocity * dt * self.sim_speed)
            e.time += dt * self.sim_speed

    
    
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
            return ((4.1357e-15)*frequency - work_func)
        
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


class ShowEField(VMobject):
    def __init__(
            self,
            max_voltage = 1.0,
            buffer_ratio= 0.8, 
            **kwargs
        ):
        super().__init__(**kwargs)
        self.voltage = 0
        self.max_voltage = max_voltage

        def scale_upd(arr):
            scaling_factor = np.clip(self.voltage, 0, self.max_voltage*buffer_ratio)/self.max_voltage
            for arrow_group in arr:
                for arrow in arrow_group:
                    arrow.scale(scaling_factor, about_edge= RIGHT)

        arrows = always_redraw(lambda: VGroup(VGroup(
            Arrow(start= (2,1.5,0), # type: ignore
                  end= (2,1.5,0) + LEFT,# type: ignore
                  buff = 0,
                  max_tip_length_to_length_ratio=0.17,
                  max_width_to_length_ratio=0.1,
                  ).shift(i*LEFT,).set_opacity(0.3).set_color(BLUE_E)# type: ignore
            for i in range(4)
            ).shift(0.5* i * DOWN)# type: ignore
            for i in range(7)

        )
        )
        arrows.add_updater(scale_upd)
        self.add(arrows)# type: ignore

