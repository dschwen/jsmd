function initRender3D( sim, container ) {
  var container;

  var renderer = null, i;
  var cam = { phi: 0, theta: 0, r: 100 };
  var intHandler = null;
  var dragging = false;

  var cw = $(container).width();
  var ch = $(container).height();

  var camera = new THREE.Camera( 30, cw/ch, 1, 10000 ); // aspect ratio 800/500!
  var scene = new THREE.Scene();
  scene.fog = new THREE.Fog( 0xffffff, 1, 10000 );

  var geometry, lines = [];

  // select renderer
  if (typeof Float32Array != "undefined") {
    try {
      renderer = new THREE.WebGLRenderer( { clearColor: 0xffffff ,clearAlpha: 1 }  );
    } catch(e) { }
  }
  if( !renderer ) {
    renderer = new THREE.CanvasRenderer( { clearColor: 0xffffff ,clearAlpha: 0 }  );
  }
  renderer.setSize( cw, ch );
  geometry = new THREE.Geometry();

  var sprite = THREE.ImageUtils.loadTexture( "ball.png" );
  var material = new THREE.ParticleBasicMaterial( { size: 1.5, map: sprite, vertexColors: false } );
  material.color.setRGB( 1.0, 0.0, 0.0 );

  // atoms
  for( i = 0; i < sim.atoms.length; ++i )
  {
    // TODO: element styling
    vector = new THREE.Vector3( sim.atoms[i].p.x-sim.ss.x/2, 
                                sim.atoms[i].p.y-sim.ss.y/2, 
                                sim.atoms[i].p.z-sim.ss.z/2 );
    geometry.vertices.push( new THREE.Vertex( vector ) );
  }

  particles = new THREE.ParticleSystem( geometry, material );
  particles.sortParticles = true;
  particles.updateMatrix();
  scene.addObject( particles );

  
  // lines
  var lm = new THREE.LineBasicMaterial( { color: 0xff0000, opacity: 0.5 , linewidth: 2} ),
      lg = [ 
        [ [1,0,0], [1,1,0] ], [ [0,1,0], [1,1,0] ], 
        [ [1,0,0], [0,0,0] ], [ [0,1,0], [0,0,0] ], 
        [ [0,0,0], [0,0,1] ], [ [1,0,0], [1,0,1] ],
        [ [0,1,0], [0,1,1] ], [ [1,1,0], [1,1,1] ],
        [ [1,0,1], [1,1,1] ], [ [0,1,1], [1,1,1] ], 
        [ [0,0,1], [0,1,1] ], [ [1,0,1], [0,0,1] ]
      ];

  for (i = 0; i < lg.length; i++) {
    geometry = new THREE.Geometry();

    vector0 = new THREE.Vector3( (lg[i][0][0]-0.5)*sim.ss.x,  
                                 (lg[i][0][1]-0.5)*sim.ss.y, 
                                 (lg[i][0][2]-0.5)*sim.ss.z );
    geometry.vertices.push( new THREE.Vertex( vector0 ) );

    vector1 = new THREE.Vector3( (lg[i][1][0]-0.5)*sim.ss.x,  
                                 (lg[i][1][1]-0.5)*sim.ss.y, 
                                 (lg[i][1][2]-0.5)*sim.ss.z );
    geometry.vertices.push( new THREE.Vertex( vector1 ) );

    lines[i] = new THREE.Line( geometry, lm ) 
    scene.addObject(lines[i]);
  }

  updateCamera();

  var light = new THREE.DirectionalLight( 0xffffff );
  light.position.x = 1;
  light.position.y = 0;
  light.position.z = 0;
  scene.addLight( light );

  $(container).empty().append( renderer.domElement );
  $(renderer.domElement).bind( 'mousemove', updateCamera  );
  $(renderer.domElement).bind( 'mousedown', function(e) { dragging = { x: e.clientX, y: e.clientY }; } );
  $(renderer.domElement).bind( 'mouseup', function() { dragging = false;  } );
  //$(renderer.domElement).bind( 'mouseout', function() { dragging = false; } );


  function updateCamera(e) {
    if( e !== undefined && dragging !== false ) {
      var dt = ( dragging.x -e.clientX ) * 0.005 ,
          df = -( dragging.y -e.clientY ) * 0.01
      camera.position.x = cam.r * Math.sin(cam.theta+dt) * Math.cos(cam.phi+df);
      camera.position.y = cam.r * Math.sin(cam.theta+dt) * Math.sin(cam.phi+df);
      camera.position.z = cam.r * Math.cos(cam.theta+dt);
    } else {
      camera.position.x = cam.r * Math.sin(cam.theta) * Math.cos(cam.phi);
      camera.position.y = cam.r * Math.sin(cam.theta) * Math.sin(cam.phi);
      camera.position.z = cam.r * Math.cos(cam.theta);
    }
  }

  function updateScene() 
  {
    var i;
    // update vertices
    for( i = 0; i < sim.atoms.length; ++i )
    {
      particles.geometry.vertices[i].position.x = sim.atoms[i].p.x - sim.ss.x/2;
      particles.geometry.vertices[i].position.y = sim.atoms[i].p.y - sim.ss.y/2;
      particles.geometry.vertices[i].position.z = sim.atoms[i].p.z - sim.ss.z/2;
    }
    particles.geometry.__dirtyVertices = true;
    renderer.render( scene, camera );
  }

  function generateSprite() {
    var canvas = document.createElement( 'canvas' );
    canvas.width = 16;
    canvas.height = 16;

    var context = canvas.getContext( '2d' );
    var gradient = context.createRadialGradient( canvas.width / 2, canvas.height / 2, 0, canvas.width / 2, canvas.height / 2, canvas.width / 2 );
    gradient.addColorStop( 0, 'rgba(255,255,255,1)' );
    gradient.addColorStop( 0.2, 'rgba(0,255,255,1)' );
    gradient.addColorStop( 0.4, 'rgba(0,0,64,1)' );
    gradient.addColorStop( 1, 'rgba(0,0,0,1)' );

    context.fillStyle = gradient;
    context.fillRect( 0, 0, canvas.width, canvas.height );

    return canvas;
  }

  return updateScene;
};
