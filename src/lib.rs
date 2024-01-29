use std::vec;
use num::complex::Complex;
use num::complex::ComplexFloat;

mod js {
    extern "C" {
        pub fn moebius_transformation(fv: Vec<f32>);
        pub fn set_distance_between_levels(distance: f32);
        pub fn set_cell_height(node: u16);
        pub fn set_children_width(width: f32);
        pub fn set_bulge(bulge: f32);
        pub fn set_base_node(node: u16);
        pub fn tree_data_sub_image(x: u32, y: u32, width: u32, height: u32, data_ptr: *const u16);
        pub fn schedule_tick();
    }
}

#[derive(Copy, Clone)]
struct Moebius {
    pub a: Complex<f32>,
    pub b: Complex<f32>,
}

impl Moebius {
    fn identity() -> Self {
        Self { 
            a: Complex::new(1.0, 0.0), 
            b: Complex::new(0.0, 0.0), 
        }
    }
}

struct TreeData {
    pub parent: u16,
    pub index: u16,
    pub prev: u16,
    pub next: u16,
    pub children: Vec<u16>,
    pub padding: Vec<u16>,
}

struct Deceleration {
    pub current: Moebius,
    pub target: Moebius,
    pub factor: f32,
}

enum PanType {
    None,
    Scroll,
    SwipeBack,
    SwipeForward,
}

static mut BASE_NODE: u16 = 0;

static mut SCROLL_AMOUNT: Vec<f32> = vec![];
static mut SCROLL_NODE: Vec<u16> = vec![];
static mut FOCUSED_SCROLL: u16 = 0;
static mut SCROLL_TARGET: f32 = 0.0;
static mut SWIPE_BACK_SCROLL_LIMIT: f32 = 0.0;

static mut PAN_START: Complex<f32> = Complex::new(0.0, 0.0);
static mut PAN_PREVIOUS: Complex<f32> = Complex::new(0.0, 0.0);
static mut PAN_PREVIOUS2: Complex<f32> = Complex::new(0.0, 0.0);
static mut PAN_CURRENT: Complex<f32> = Complex::new(0.0, 0.0);

static mut DURATION_BETWEEN_PAN_PREVIOUS_AND_CURRENT: f32 = 0.0;
static mut DURATION_BETWEEN_PAN_PREVIOUS2_AND_CURRENT: f32 = 0.0;
static mut DURATION_SINCE_PAN: f32 = 0.0;

static mut MOUSE_DOWN_X: f32 = 0.0;
static mut MOUSE_DOWN_Y: f32 = 0.0;

static mut WILL_PAN: bool = false;
static mut INTERACTION_IS_TOUCH: bool = false;
static mut INTERACTION_STOPPED_SCROLLING: bool = false;

static mut PAN_TYPE: PanType = PanType::None;

pub struct App {
    viewport_width: f32,
    viewport_height: f32,
    tree_data: Vec<TreeData>,

    distance_between_levels: f32,
    cell_height: f32,
    children_width: f32,
    buldge: f32,

    number_of_nodes: u16,
    focused_note: u16,
}

impl App {
    fn moebius_equal(m1: Moebius, m2: Moebius) -> bool {
        return m1.a == m2.a && m1.b == m2.b;
    }

    fn compose(m1: Moebius, m2: Moebius) -> Moebius {
        let c = m2.a + m1.b * m2.b.conj();
        
        Moebius { 
            a: m1.a * c / (m2.a * m2.b * m1.b + 1.0),
            b: (m2.a * m2.b + m1.b) / c,
        }
    }

    fn inverse(m: Moebius) -> Moebius {
        Moebius { 
            a: m.a.conj(),
            b: -m.b * m.a ,
        }
    }
    
    fn apply(m: Moebius, z: Complex<f32>) -> Complex<f32> {
        m.a * (z + m.b) / (m.b.conj() * z + 1.0)
    }

    fn origin_translation(z: Complex<f32>) -> Moebius {
        Moebius {
            a: Complex::new(1.0, 0.0),
            b: z,
        }
    }

    fn mpow(m: Moebius, k: f32) -> Moebius {
        if (m.a - 1.0).abs() < 1e-6 {
            let len = m.b.abs();
    
            if len < 1e-6 {
                return Moebius::identity()
            }
    
            let xk = f32::powf(1.0 - len, k);
            let yk = f32::powf(1.0 + len, k);
    
            Moebius {
                a: Complex { re: 1.0, im: 0.0 },
                b: m.b / len * (yk - xk) / (xk + yk),
            }
        } else {
            let d = Complex::sqrt(4.0 * m.a * m.b * m.b.conj() + (m.a - 1.0) * (m.a - 1.0));
            let a = d - m.a + 1.0;
            let b = d + m.a - 1.0;
    
            let xk = (-d + m.a + 1.0).powf(k);
            let yk = (d + m.a + 1.0).powf(k);
    
            Moebius {
                a: (a * xk + b * yk) / (b * xk + a * yk),
                b: m.a * m.b * 2.0 * (yk - xk) / (a * xk + b * yk),
            }
        }
    }

    fn disk_to_band(&self, disk: Complex<f32>) -> Complex<f32> {
        (disk * Complex::i()).atanh() / (Complex::i() * 0.785) * self.viewport_width
    }

    fn band_to_disk(&self, band: Complex<f32>) -> Complex<f32> {
        (band / self.viewport_width * (Complex::i() * 0.785)).tanh() / Complex::i()
    }

    fn row_in_band(&self, band: Complex<f32>) -> i32 {
        ((band.im + self.viewport_height * 0.5) / self.cell_height).floor() as i32
    }

    fn scroll_in_node(&self, scroll: f32, node: &mut u16) -> f32 {
        let mut offset: f32 = 0.0;

        for _ in 1..32 {
            let data = &self.tree_data[*node as usize];

            if scroll + offset < 0.0 {
                if data.prev == 0 {
                    break
                }

                *node = data.prev;
            } else if scroll + offset >= self.cell_height * 4.0 {
                if data.prev == 0 {
                    break
                }

                *node = data.next;
                offset -= self.cell_height * 4.0;
            } else {
                break
            }
        }

        offset
    }

    fn move_from_parent(&self, row: f32, row_in_parent: f32) -> Moebius {
        App::compose(
            App::origin_translation(
                self.band_to_disk(Complex::i() * (self.cell_height * (row_in_parent + 0.5) - self.viewport_height * 0.5)),
            ),
            App::compose(
                App::origin_translation(Complex::new(self.distance_between_levels, 0.0)), 
                App::origin_translation(self.band_to_disk(Complex::i() * self.cell_height * row)),
            )
        )
    }

    // Assumes that parent nodes always have a lower node number than their children
    // and that earlier nodes in a list always have a lower node number than later
    // nodes.
    fn move_between_nodes(&self, mut from: u16, mut to: u16) -> Moebius {
        let mut moebius_from = Moebius::identity();
        let mut moebius_to = Moebius::identity();

        while from != to {
            if self.tree_data[to as usize].parent == self.tree_data[from as usize].parent {
                moebius_to = App::compose(
                    App::origin_translation(
                        self.band_to_disk(
                            Complex::i() * self.cell_height * (self.tree_data[to as usize].index as f32)
                        )
                    ), 
                    moebius_to,
                );
                to = from;
            } else if from < to {
                moebius_to = App::compose(
                    self.move_from_parent(
                        self.tree_data[to as usize].index as f32, 
                        (self.tree_data[to as usize].parent as f32).rem_euclid(4.0),
                    ),
                    moebius_to,
                );
                to = self.tree_data[to as usize].parent / 4;
            } else {
                moebius_from = App::compose(
                    self.move_from_parent(
                        self.tree_data[from as usize].index as f32, 
                        (self.tree_data[from as usize].parent as f32).rem_euclid(4.0),
                    ), 
                    moebius_from,
                );
                from = self.tree_data[from as usize].parent / 4;
            }
        }

        App::compose(App::inverse(moebius_from), moebius_to)
    }

    fn domain_tranformation_step(&self, z: Complex<f32>, node: &mut u16) -> Moebius {
        let data = &self.tree_data[*node as usize];
        let parent = data.parent / 4;
        let mut band = self.disk_to_band(z);
        let parent_band = self.disk_to_band(
            App::apply(
                App::origin_translation(self.distance_between_levels.into()), 
                self.band_to_disk(band + Complex::i() * self.cell_height * (data.index as f32))
            )
        ) + Complex::i() * 0.5 * self.cell_height;

        if parent_band.re < self.viewport_width - self.children_width || (parent_band.im - 0.5 * self.cell_height).abs() > 0.5 * self.cell_height {
            if parent == 0 {
                return Moebius::identity()
            }

            *node = parent;

            return self.move_from_parent(data.index as f32, (data.parent % 4) as f32)
        } else {
            let scroll_offset = self.scroll_in_node(band.im + self.viewport_height * 0.5, node);

            if scroll_offset != 0.0 {
                return App::origin_translation(self.band_to_disk(Complex::i() * scroll_offset))
            }

            let row = self.row_in_band(band);
            band -= Complex::i() * ((row as f32) * self.cell_height - self.viewport_height * 0.5);

            let mut next: u16 = 0;

            if row >= 0 && row < 4 {
                next = data.children[row as usize];
            }

            if band.re < self.viewport_width - self.children_width || next == 0 {
                return Moebius::identity()
            }

            *node = next;

            App::inverse(self.move_from_parent(0 as f32, row as f32))
        }
    }

    fn domain_tranformation(&self, z: Complex<f32>, node: &mut u16) -> Moebius {
        let mut m = Moebius::identity();

        for _ in 1..32 {
            let d = self.domain_tranformation_step(App::apply(m, z), node);

            if App::moebius_equal(d, Moebius::identity()) {
                break
            }

            m = App::compose(d, m);
        }

        m
    }

    fn rubberband(&self, scroll_amount: f32, constrained: f32, factor: f32) -> f32 {
        let delta = constrained - scroll_amount;

        if delta == 0.0 {
            return scroll_amount
        }

        constrained - f32::copysign(
            (1.0 - (1.0 / (delta.abs() * factor / self.viewport_height) + 1.0)) * self.viewport_height, 
            delta
        )
    }

    fn unrubberband(&self, scroll_amount: f32, constrained: f32, factor: f32) -> f32 {
        let delta = constrained - scroll_amount;

        if delta == 0.0 {
            return scroll_amount
        }

        constrained - f32::copysign(
            self.viewport_height / factor * (1.0 / (1.0 - delta.abs() / self.viewport_height) - 1.0), 
            delta
        )
    }

    fn falloff(x: f32) -> f32 {
        // x.tanh()
        1.0 - 1.0 / x.exp()
    }

    fn pan_location_for_xy(&self, mut x: f32, y: f32, has_parent: bool) -> Complex<f32> {
        if !has_parent && x > 0.0 {
            x = self.rubberband(x, 0.0, 0.4);
        }

        let ramp: f32 = 0.25;

        if x / self.viewport_width >= 1.0 - ramp {
            x = self.viewport_width * (1.0 + ramp * (App::falloff((x / self.viewport_width - 1.0 + ramp) / ramp) - 1.0));
        }

        x *= 0.9;

        self.band_to_disk(x + Complex::i() * y)
    }

    fn moebius_for_pan(&self, start: Complex<f32>, current: Complex<f32>) -> Moebius {
        let m = App::origin_translation(-start);

        let delta = App::apply(
            App::inverse(m), 
            -current
        );

        App::compose(
            App::inverse(m), 
            App::compose(App::origin_translation(delta), m),
        )
    }

    fn set_moebius_transformation(m: Moebius) {
        unsafe {
            js::moebius_transformation(vec![m.a.re(), m.a.im(), m.b.re(), m.b.im()])
        }
    }

    fn deceleration_active(d: &mut Deceleration) -> bool {
        !App::moebius_equal(d.target, Moebius::identity()) || 
        !App::moebius_equal(d.current, Moebius::identity())
    }

    fn decelerate(d: &mut Deceleration, dt: f32) -> bool {
        if !App::deceleration_active(d) {
            return false
        }

        let k = d.factor.powf(1000.0 * dt);

        d.current = App::compose(
            App::mpow(App::compose(d.current, App::inverse(d.target)), k), 
            d.target,
        );

        d.target = App::mpow(d.target, k);

        d.current.a /= d.current.a.abs();
        d.target.a /= d.target.a.abs();

        if d.current.b.abs() > 1.0 {
            d.current.b /= d.current.b.abs();
        }

        if d.target.b.abs() > 1.0 {
            d.target.b /= d.target.b.abs();
        }

        if App::deceleration_active(d) {
            unsafe { 
                js::schedule_tick(); 
            }
        }

        true
    }

    fn deceleration_offset_for_velocity(dm: Moebius, dt: f32, factor: f32) -> Moebius {
        if dt <= 0.0 || App::moebius_equal(dm, Moebius::identity()) {
            return Moebius::identity();
        }
        let k = f32::powf(factor, dt * 1000.0);
        
        return App::mpow(dm, 1.0 / (1.0 - k));
    }
    
    fn constrained_scroll_amount(&self, scroll_amount: f32, scroll_node: u16) -> f32 {
        let mut scroll_amount: f32 = scroll_amount;
        let mut scroll_node = scroll_node;
    
        let bottom_offset = self.scroll_in_node(scroll_amount + self.viewport_height, &mut scroll_node);
    
        if scroll_amount + bottom_offset + self.viewport_height > 4.0 * self.cell_height {
            scroll_amount = 4.0 * self.cell_height - bottom_offset - self.viewport_height;
        }

        let top_offset = self.scroll_in_node(scroll_amount, &mut scroll_node);
    
        if scroll_amount + top_offset < 0.0 {
            scroll_amount = -top_offset;
        }
    
        return scroll_amount;
    }
}

#[no_mangle]
pub extern "C" fn start() {
    unsafe {
        //let app = App::new();

        //js::moebius_transformation();
        //js::clearToBlue();
    }
}